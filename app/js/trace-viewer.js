const TraceViewer = {
    init(manifest) {
        this.container = document.getElementById('trace-container');
        this.toolLog = manifest.tool_log || [];
        this.render();
    },

    render() {
        this.container.innerHTML = '';
        if (this.toolLog.length === 0) {
            this.container.innerHTML = '<div class="empty">No tool calls recorded.</div>';
            return;
        }

        const list = document.createElement('div');
        list.className = 'trace-list';

        this.toolLog.forEach((entry, idx) => {
            const item = document.createElement('div');
            item.className = `trace-item status-${entry.status}`;
            item.dataset.idx = idx;
            
            // Format arguments nicely
            let argsHtml = '';
            if (entry.arguments) {
                if (entry.name === 'run_python' && entry.arguments.code) {
                    argsHtml = `<div class="trace-args"><strong>Code:</strong><pre class="code-block">${this.escapeHtml(entry.arguments.code)}</pre></div>`;
                } else {
                    // Filter out large content strings for other tools
                    const cleanArgs = {...entry.arguments};
                    if (cleanArgs.content && cleanArgs.content.length > 500) {
                        cleanArgs.content = cleanArgs.content.substring(0, 500) + '... [truncated]';
                    }
                    argsHtml = `<div class="trace-args"><strong>Params:</strong><pre class="json-block">${this.escapeHtml(JSON.stringify(cleanArgs, null, 2))}</pre></div>`;
                }
            }
            
            const timestamp = entry.timestamp ? new Date(entry.timestamp).toLocaleTimeString([], { hour: '2-digit', minute: '2-digit', second: '2-digit' }) : '';

            item.innerHTML = `
                <div class="trace-header">
                    <span class="idx">#${idx}</span>
                    <span class="tool-name">${entry.name}</span>
                    ${timestamp ? `<span class="timestamp-pill">${timestamp}</span>` : ''}
                    <span class="duration">${entry.duration_ms}ms</span>
                    <span class="status-badge">${entry.status}</span>
                </div>
                <div class="trace-preview">${this.escapeHtml(entry.preview)}</div>
                <div class="trace-details hidden">
                    ${entry.reasoning ? `<div class="trace-reasoning"><strong>Reasoning:</strong><p>${this.escapeHtml(entry.reasoning)}</p></div>` : ''}
                    ${argsHtml}
                    ${entry.result ? `<div class="trace-result"><strong>Full Result:</strong><pre class="json-block">${this.escapeHtml(entry.result)}</pre></div>` : ''}
                    ${entry.tokens ? `<div class="trace-tokens"><strong>Tokens:</strong> <span class="badge">${entry.tokens.total} total</span> <small>(${entry.tokens.prompt}p + ${entry.tokens.completion}c)</small></div>` : ''}
                    
                    <div class="step-feedback-box">
                        <strong>Feedback for Step #${idx} (Updates Agent Prompt):</strong>
                        <textarea 
                            class="feedback-input" 
                            id="feedback-step-${idx}"
                            placeholder="Describe what the agent should have done differently. Your feedback will be used to improve the agent's prompt for future runs."
                        ></textarea>
                        <div class="step-feedback-actions">
                            <button class="small-btn primary" onclick="TraceViewer.submitStepFeedback(${idx})">Submit Feedback</button>
                        </div>
                    </div>
                </div>
            `;

            item.onclick = (e) => {
                // Prevent bubbling if clicking feedback elements
                if (e.target.closest('.step-feedback-box')) return;
                
                this.toggleDetails(item, entry);
                
                // Highlight the trace item itself
                this.highlight(idx);
                
                // Trigger DAG highlighting if available
                if (window.DagRenderer) {
                    DagRenderer.highlightForStep(idx);
                }
                
                // Trigger Rubric highlighting
                if (window.RubricViewer) {
                    RubricViewer.highlightForStep(idx);
                }
            };

            list.appendChild(item);
        });

        this.container.appendChild(list);
    },

    toggleDetails(item, entry) {
        const details = item.querySelector('.trace-details');
        details.classList.toggle('hidden');
        item.classList.toggle('expanded');
    },

    escapeHtml(text) {
        const div = document.createElement('div');
        div.textContent = text;
        return div.innerHTML;
    },

    highlight(idx) {
        const items = this.container.querySelectorAll('.trace-item');
        items.forEach(item => item.classList.remove('highlighted'));
        const target = this.container.querySelector(`.trace-item[data-idx="${idx}"]`);
        if (target) {
            target.classList.add('highlighted');
            target.scrollIntoView({ behavior: 'smooth', block: 'nearest' });
        }
    },

    submitStepFeedback(idx) {
        const input = this.container.querySelector(`#feedback-step-${idx}`);
        if (!input) return;
        
        const suggestion = input.value;
        if (suggestion && suggestion.trim() !== '') {
            FeedbackSystem.submitDirect('action', `step:${idx}`, suggestion);
            
            // Visual feedback
            const btn = input.parentElement.querySelector('button');
            const originalText = btn.textContent;
            btn.textContent = 'Saved!';
            btn.classList.add('success');
            setTimeout(() => {
                btn.textContent = originalText;
                btn.classList.remove('success');
            }, 2000);
        } else {
            alert('Please enter some feedback first.');
        }
    }
};

