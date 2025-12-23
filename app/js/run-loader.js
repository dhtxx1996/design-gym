const RunLoader = {
    async listTasks() {
        console.log('Listing tasks from index...');
        try {
            const response = await fetch('data/index.json');
            if (!response.ok) throw new Error('Index not found');
            const data = await response.json();
            this.cachedIndex = data;
            return Object.keys(data);
        } catch (e) {
            console.warn('Could not load index.json, using fallback', e);
            return ['ph_sensitive_design']; 
        }
    },

    async listRuns(task) {
        console.log(`Listing runs for task: ${task}`);
        try {
            if (!this.cachedIndex) {
                await this.listTasks();
            }
            return this.cachedIndex[task] || [];
        } catch (e) {
            console.error(`Error listing runs for ${task}:`, e);
            return [];
        }
    },

    async loadManifest(task, runId) {
        const path = `data/${task}/${runId}/manifest.json`;
        console.log(`Fetching manifest from: ${path}`);
        try {
            const response = await fetch(path);
            if (!response.ok) {
                console.error(`Fetch failed for ${path}: ${response.status} ${response.statusText}`);
                return null;
            }
            const data = await response.json();
            console.log(`Successfully loaded manifest for ${runId}`);
            return data;
        } catch (e) {
            console.error(`Error loading manifest from ${path}:`, e);
            return null;
        }
    }
};

