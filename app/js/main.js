// App initialization and global state
document.addEventListener('DOMContentLoaded', () => {
    console.log('App initialization started');
    const path = window.location.pathname;
    const urlParams = new URLSearchParams(window.location.search);
    const runId = urlParams.get('runId');
    const task = urlParams.get('task');

    if (path.includes('run-details.html') && runId && task) {
        console.log('Routing to Run Details');
        initRunDetails(task, runId);
    } else {
        console.log('Routing to Dashboard');
        initDashboard();
    }
});

function initDashboard() {
    console.log('initDashboard called');
    if (typeof Dashboard !== 'undefined') {
        Dashboard.init();
    } else {
        console.error('Dashboard object not found! check if dashboard.js is loaded.');
        // Retry once after a short delay
        setTimeout(() => {
            if (typeof Dashboard !== 'undefined') {
                console.log('Dashboard found on retry');
                Dashboard.init();
            } else {
                const errDisplay = document.getElementById('error-display');
                if (errDisplay) {
                    errDisplay.textContent = 'Critical Error: Dashboard scripts failed to load.';
                    errDisplay.classList.remove('hidden');
                }
            }
        }, 500);
    }
}

async function initRunDetails(task, runId) {
    const manifest = await RunLoader.loadManifest(task, runId);
    if (!manifest) {
        alert('Could not load manifest for ' + runId);
        return;
    }

    document.getElementById('run-id-title').textContent = `Run: ${runId}`;
    
    // Initialize components
    if (typeof TraceViewer !== 'undefined') TraceViewer.init(manifest);
    if (typeof DagRenderer !== 'undefined') DagRenderer.init(manifest);
    if (typeof RubricViewer !== 'undefined') RubricViewer.init(manifest);
    if (typeof FeedbackSystem !== 'undefined') FeedbackSystem.init(manifest);
    if (typeof OptimizationEngine !== 'undefined') OptimizationEngine.init(manifest);
    
    // Toggle logic
    const serialBtn = document.getElementById('toggle-serial');
    const dagBtn = document.getElementById('toggle-dag');
    const serialView = document.getElementById('serial-view');
    const dagView = document.getElementById('dag-view');

    serialBtn.onclick = () => {
        serialBtn.classList.add('active');
        dagBtn.classList.remove('active');
        serialView.classList.remove('hidden');
        dagView.classList.add('hidden');
    };

    dagBtn.onclick = () => {
        dagBtn.classList.add('active');
        serialBtn.classList.remove('active');
        dagView.classList.remove('hidden');
        serialView.classList.add('hidden');
        DagRenderer.render(); // Ensure render on first show
    };

    // Initialize Resizer
    initResizer();

    renderMetadata(manifest);
}

function initResizer() {
    const layout = document.querySelector('.details-layout');
    const resizer = document.getElementById('drag-handle');
    if (!layout || !resizer) return;

    let isResizing = false;

    resizer.addEventListener('mousedown', (e) => {
        isResizing = true;
        resizer.classList.add('active');
        document.body.style.cursor = 'col-resize';
        // Prevent text selection during drag
        document.body.style.userSelect = 'none';
    });

    document.addEventListener('mousemove', (e) => {
        if (!isResizing) return;
        
        // Calculate new width for the right panel
        // The right panel is on the right side, so we measure from the right edge of the window
        const newWidth = window.innerWidth - e.clientX;
        
        // Constrain width (min 300px, max 800px)
        if (newWidth > 300 && newWidth < 800) {
            // Update grid template columns
            // 280px left panel, 1fr center, 4px resizer, newWidth right panel
            layout.style.gridTemplateColumns = `280px 1fr 4px ${newWidth}px`;
        }
    });

    document.addEventListener('mouseup', () => {
        if (isResizing) {
            isResizing = false;
            resizer.classList.remove('active');
            document.body.style.cursor = '';
            document.body.style.userSelect = '';
        }
    });
}

function renderMetadata(manifest) {
    const container = document.querySelector('#run-metadata .meta-content');
    const exec = manifest.execution || {};
    
    container.innerHTML = `
        <div class="meta-item"><strong>Task:</strong> ${manifest.task}</div>
        <div class="meta-item"><strong>Model:</strong> ${exec.agent_model}</div>
        <div class="meta-item"><strong>Duration:</strong> ${exec.duration_seconds}s</div>
        <div class="meta-item"><strong>Tokens:</strong> ${exec.tokens?.total || 0}</div>
        <div class="meta-item"><strong>Completed:</strong> ${exec.task_completed ? '✅' : '❌'}</div>
    `;

    const sourceList = document.querySelector('.source-list');
    if (manifest.prompts) {
        Object.entries(manifest.prompts).forEach(([key, val]) => {
            if (key.endsWith('_source')) {
                const li = document.createElement('li');
                li.className = 'source-item';
                li.textContent = `${key.replace('_source', '')}: ${val}`;
                sourceList.appendChild(li);
            }
        });
    }
}

