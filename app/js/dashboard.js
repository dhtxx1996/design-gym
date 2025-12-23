const Dashboard = {
    async init() {
        console.log('Dashboard initialization...');
        this.runGrid = document.getElementById('run-grid');
        this.taskFilter = document.getElementById('task-filter');
        this.searchInput = document.getElementById('search-input');
        
        if (!this.runGrid || !this.taskFilter || !this.searchInput) {
            console.error('Critical UI elements missing from dashboard');
            return;
        }
        
        try {
            const tasks = await RunLoader.listTasks();
            this.renderTaskFilter(tasks);
            
            await this.loadAndRenderAll();
            
            this.taskFilter.onchange = () => this.loadAndRenderAll();
            this.searchInput.oninput = () => this.filterRuns();
            
            this.startPolling();
            console.log('Dashboard initialized successfully');
        } catch (e) {
            console.error('Failed to initialize dashboard:', e);
            if (this.runGrid) {
                this.runGrid.innerHTML = `<div class="error">Failed to initialize dashboard: ${e.message}</div>`;
            }
        }
    },

    renderTaskFilter(tasks) {
        tasks.forEach(task => {
            const opt = document.createElement('option');
            opt.value = task;
            opt.textContent = task;
            this.taskFilter.appendChild(opt);
        });
    },

    async loadAndRenderAll() {
        console.log('Loading all runs...');
        this.runGrid.innerHTML = '<div class="loading">Loading runs...</div>';
        
        try {
            const selectedTask = this.taskFilter.value;
            const tasks = selectedTask ? [selectedTask] : await RunLoader.listTasks();
            console.log('Tasks to load:', tasks);
            
            let allRuns = [];
            for (const task of tasks) {
                const runs = await RunLoader.listRuns(task);
                console.log(`Runs for task ${task}:`, runs);
                for (const runId of runs) {
                    const manifest = await RunLoader.loadManifest(task, runId);
                    if (manifest) {
                        allRuns.push(manifest);
                    } else {
                        console.warn(`Failed to load manifest for ${task}/${runId}`);
                    }
                }
            }
            
            console.log('Total runs loaded:', allRuns.length);
            this.allManifests = allRuns;
            
            if (allRuns.length === 0) {
                this.runGrid.innerHTML = '<div class="empty">No runs found. Check if outputs directory contains manifest.json files and if data/ symlink is correct.</div>';
            } else {
                this.renderRuns(allRuns);
            }
            this.renderStats(allRuns);
        } catch (e) {
            console.error('Error in loadAndRenderAll:', e);
            this.runGrid.innerHTML = `<div class="error">Error loading runs: ${e.message}</div>`;
        }
    },

    renderRuns(manifests) {
        this.runGrid.innerHTML = '';
        manifests.forEach(m => {
            const card = document.createElement('div');
            card.className = 'run-card';
            const exec = m.execution || {};
            const statusClass = exec.task_completed ? 'status-completed' : 'status-failed';
            
            card.innerHTML = `
                <span class="status ${statusClass}">${exec.task_completed ? 'Completed' : 'Failed'}</span>
                <h3>${m.run_id}</h3>
                <div class="task-name">${m.task}</div>
                <div class="metrics">
                    <span class="badge">Model: ${exec.agent_model}</span>
                    <span class="badge">Duration: ${exec.duration_seconds}s</span>
                </div>
            `;
            card.onclick = () => {
                window.location.href = `run-details.html?task=${m.task}&runId=${m.run_id}`;
            };
            this.runGrid.appendChild(card);
        });
    },

    renderStats(manifests) {
        const statsPanel = document.getElementById('stats-panel');
        if (!statsPanel) return;

        const total = manifests.length;
        const completed = manifests.filter(m => m.execution?.task_completed).length;
        const avgDuration = manifests.reduce((acc, m) => acc + (m.execution?.duration_seconds || 0), 0) / (total || 1);
        
        statsPanel.innerHTML = `
            <div class="stat-card">
                <span class="stat-label">Total Runs</span>
                <span class="stat-value">${total}</span>
            </div>
            <div class="stat-card">
                <span class="stat-label">Success Rate</span>
                <span class="stat-value">${Math.round((completed / (total || 1)) * 100)}%</span>
            </div>
            <div class="stat-card">
                <span class="stat-label">Avg Duration</span>
                <span class="stat-value">${Math.round(avgDuration)}s</span>
            </div>
            <div class="stat-card">
                <span class="stat-label">Tasks</span>
                <span class="stat-value">${new Set(manifests.map(m => m.task)).size}</span>
            </div>
        `;
    },

    filterRuns() {
        const term = this.searchInput.value.toLowerCase();
        const filtered = this.allManifests.filter(m => 
            m.run_id.toLowerCase().includes(term) || 
            m.task.toLowerCase().includes(term)
        );
        this.renderRuns(filtered);
    },

    startPolling() {
        console.log('Starting live polling...');
        setInterval(async () => {
            // Only refresh if search and filter are empty to avoid interrupting user
            if (this.searchInput.value === '' && this.taskFilter.value === '') {
                await this.refreshSilently();
            }
        }, 10000); // Poll every 10 seconds
    },

    async refreshSilently() {
        const tasks = await RunLoader.listTasks();
        let allRuns = [];
        for (const task of tasks) {
            const runs = await RunLoader.listRuns(task);
            for (const runId of runs) {
                const manifest = await RunLoader.loadManifest(task, runId);
                if (manifest) allRuns.push(manifest);
            }
        }
        
        // Only re-render if run count or content changed
        const currentIds = (this.allManifests || []).map(m => m.run_id).join(',');
        const newIds = allRuns.map(m => m.run_id).join(',');
        
        if (currentIds !== newIds) {
            console.log('New runs detected! Refreshing dashboard.');
            this.allManifests = allRuns;
            this.renderRuns(allRuns);
            this.renderStats(allRuns);
        }
    }
};

