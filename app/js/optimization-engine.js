const OptimizationEngine = {
    init(manifest) {
        this.container = document.getElementById('suggestions-container');
        this.manifest = manifest;
        this.suggestions = [];
        this.generateInitialSuggestions();
    },

    generateInitialSuggestions() {
        // Automatically generate some suggestions based on the manifest (e.g. failures)
        const exec = this.manifest.execution || {};
        const toolCalls = this.manifest.tool_log || [];
        const failures = toolCalls.filter(t => t.status === 'error');

        if (failures.length > 0) {
            this.addSuggestion({
                type: 'prompt_tuning',
                priority: 'high',
                title: 'Address tool usage errors',
                description: `Agent failed ${failures.length} tool calls. Consider adding clearer constraints in the system prompt for tools: ${[...new Set(failures.map(f => f.name))].join(', ')}.`
            });
        }

        const score = this.manifest.evaluations?.[0]?.score?.normalized;
        if (score !== undefined && score < 0.8) {
            this.addSuggestion({
                type: 'rubric_refinement',
                priority: 'medium',
                title: 'Review low-scoring criteria',
                description: 'The overall score is below 80%. Inspect category breakdowns to refine agent instructions or the rubric itself for low-performing areas.'
            });
        }

        // Add specific rubric suggestions if they have low scores
        const latestEval = this.manifest.evaluations?.[this.manifest.evaluations.length - 1];
        if (latestEval && latestEval.categories) {
            Object.entries(latestEval.categories).forEach(([key, cat]) => {
                const ratio = cat.mean / (cat.max || 1);
                if (ratio < 0.5) {
                    this.addSuggestion({
                        type: 'rubric_refinement',
                        priority: 'high',
                        title: `Critical Failure: ${key}`,
                        description: `Criterion "${key}" scored poorly (${cat.mean}/${cat.max}). The reasoning says: "${Array.isArray(cat.reasoning) ? cat.reasoning[0] : cat.reasoning}". Consider if the rubric is too strict or the agent prompt lacks instructions for this.`
                    });
                }
            });
        }

        this.render();
    },

    onNewFeedback(feedback) {
        // Turn feedback into an optimization suggestion
        this.addSuggestion({
            type: feedback.category,
            priority: 'high',
            title: `Feedback-driven: ${feedback.category}`,
            description: feedback.suggestion,
            context: feedback.target
        });
        this.render();
    },

    addSuggestion(suggestion) {
        this.suggestions.unshift(suggestion); // Newest first
    },

    render() {
        this.container.innerHTML = '';
        if (this.suggestions.length === 0) {
            this.container.innerHTML = '<div class="empty">No suggestions yet. Provide feedback to generate some!</div>';
            return;
        }

        this.suggestions.forEach(s => {
            const card = document.createElement('div');
            card.className = `suggestion-card priority-${s.priority}`;
            card.innerHTML = `
                <div class="suggestion-header">
                    <span class="badge">${s.type}</span>
                    <span class="priority-label">${s.priority} priority</span>
                </div>
                <h5>${s.title}</h5>
                <p>${s.description}</p>
                ${s.context ? `<div class="context">Target: ${s.context}</div>` : ''}
                <button class="apply-btn" onclick="alert('In a real app, this would apply the optimization to your source code or prompt files.')">Draft Optimization</button>
            `;
            this.container.appendChild(card);
        });
    }
};

