const RubricViewer = {
    init(manifest) {
        this.container = document.getElementById('rubric-container');
        this.manifest = manifest;
        this.render();
    },

    getStepIndex(stepRef) {
        if (typeof stepRef === 'number') return stepRef;
        if (typeof stepRef === 'string') {
            const match = stepRef.match(/_(\d+)$/);
            if (match) return parseInt(match[1]);
            if (!isNaN(stepRef)) return Number(stepRef);
        }
        return -1;
    },

    render() {
        const rubric = this.manifest.rubric || {};
        const evaluations = this.manifest.evaluations || [];
        const latestEval = evaluations[evaluations.length - 1] || {};
        const categories = latestEval.categories || {};
        const criteria = rubric.criteria || [];

        if (criteria.length === 0) {
            this.container.innerHTML = '<div class="empty">No rubric data available.</div>';
            return;
        }

        let html = '<div class="rubric-list">';
        
        criteria.forEach(criterion => {
            const catData = categories[criterion.key] || {};
            const score = catData.mean !== undefined ? catData.mean : '?';
            const max = criterion.max || '?';
            const reasoning = Array.isArray(catData.reasoning) ? catData.reasoning[0] : catData.reasoning;
            
            const scoreClass = this.getScoreClass(score, max);

            html += `
                <div class="rubric-item" data-key="${criterion.key}" 
                     onmouseenter="RubricViewer.highlightRelated('${criterion.key}')"
                     onmouseleave="RubricViewer.clearHighlights()">
                    <div class="rubric-header">
                        <span class="criterion-name">${criterion.name}</span>
                        <span class="criterion-score ${scoreClass}">${score}/${max}</span>
                    </div>
                    ${reasoning ? `<div class="criterion-reasoning">"${reasoning}"</div>` : ''}
                    <div class="rubric-links">
                        ${this.renderLinks(catData)}
                    </div>
                    <div class="rubric-actions">
                        <textarea 
                            class="rubric-feedback-area" 
                            placeholder="Add feedback to improve grading for this criterion..."
                            onblur="FeedbackSystem.submitDirect('rubric', 'criterion:${criterion.key}', this.value)"
                        ></textarea>
                    </div>
                </div>
            `;
        });

        html += '</div>';
        
        if (rubric.raw_content) {
            html += `
                <details class="raw-rubric">
                    <summary>View Raw Rubric Text</summary>
                    <pre>${this.escapeHtml(rubric.raw_content)}</pre>
                </details>
            `;
        }

        this.container.innerHTML = html;
    },

    getScoreClass(score, max) {
        if (score === '?' || max === '?') return '';
        const ratio = score / max;
        if (ratio >= 0.8) return 'score-high';
        if (ratio >= 0.5) return 'score-med';
        return 'score-low';
    },

    renderLinks(catData) {
        let html = '';
        if (catData.related_files && catData.related_files.length > 0) {
            html += '<div class="related-links">Files: ' + 
                catData.related_files.map(f => `<span class="link-badge">${f}</span>`).join(' ') + 
                '</div>';
        }
        if (catData.related_steps && catData.related_steps.length > 0) {
            html += '<div class="related-links">Steps: ' + 
                catData.related_steps.map(s => {
                    const idx = this.getStepIndex(s);
                    return `<span class="link-badge" onclick="TraceViewer.highlight(${idx})">#${idx}</span>`;
                }).join(' ') + 
                '</div>';
        }
        return html;
    },

    highlightRelated(key) {
        const evaluations = this.manifest.evaluations || [];
        const latestEval = evaluations[evaluations.length - 1] || {};
        const catData = latestEval.categories?.[key] || {};

        if (catData.related_files) {
            catData.related_files.forEach(file => {
                const nodeId = `file:${file}`;
                if (window.DagRenderer) {
                    window.DagRenderer.highlightNode(nodeId);
                }
            });
        }

        if (catData.related_steps) {
            catData.related_steps.forEach(stepRef => {
                const stepIdx = this.getStepIndex(stepRef);
                if (window.TraceViewer) {
                    window.TraceViewer.highlight(stepIdx);
                }
                if (window.DagRenderer) {
                    const nodes = this.manifest.dag?.nodes || [];
                    const pattern = `_${stepIdx}`;
                    const node = nodes.find(n => n.id.endsWith(pattern));
                    if (node) window.DagRenderer.highlightNode(node.id);
                }
            });
        }
    },

    highlightForStep(stepIdx) {
        console.log('Highlighting rubric for step:', stepIdx);
        this.clearHighlights();
        
        const evaluations = this.manifest.evaluations || [];
        const latestEval = evaluations[evaluations.length - 1] || {};
        const categories = latestEval.categories || {};
        
        let found = false;
        Object.entries(categories).forEach(([key, data]) => {
            // Robust matching: extract index from string references (e.g. "run_python_1" -> 1)
            const steps = (data.related_steps || []).map(s => this.getStepIndex(s));
            if (steps.includes(Number(stepIdx))) {
                const el = this.container.querySelector(`.rubric-item[data-key="${key}"]`);
                if (el) {
                    console.log('Match found for rubric item:', key);
                    el.classList.add('highlighted');
                    // Scroll the rubric item into view within the right panel
                    el.scrollIntoView({ behavior: 'smooth', block: 'center' });
                    found = true;
                }
            }
        });

        // Also scroll the whole rubric section into view if any match found
        if (found) {
            const rubricSection = document.getElementById('rubric-section');
            if (rubricSection) {
                rubricSection.scrollIntoView({ behavior: 'smooth', block: 'start' });
            }
        } else {
            console.log('No rubric matches for step:', stepIdx);
        }
    },

    highlightForFile(fileName) {
        console.log('Highlighting rubric for file:', fileName);
        this.clearHighlights();
        const evaluations = this.manifest.evaluations || [];
        const latestEval = evaluations[evaluations.length - 1] || {};
        const categories = latestEval.categories || {};
        
        Object.entries(categories).forEach(([key, data]) => {
            if (data.related_files && data.related_files.includes(fileName)) {
                const el = this.container.querySelector(`.rubric-item[data-key="${key}"]`);
                if (el) {
                    console.log('Match found for rubric item:', key);
                    el.classList.add('highlighted');
                    el.scrollIntoView({ behavior: 'smooth', block: 'nearest' });
                }
            }
        });
    },

    clearHighlights() {
        if (window.DagRenderer) window.DagRenderer.clearHighlights();
        const items = this.container.querySelectorAll('.rubric-item');
        items.forEach(item => item.classList.remove('highlighted'));
    },

    escapeHtml(text) {
        const div = document.createElement('div');
        div.textContent = text;
        return div.innerHTML;
    }
};

