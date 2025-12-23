const FeedbackSystem = {
    init(manifest) {
        this.container = document.getElementById('feedback-form-container');
        this.manifest = manifest;
        this.feedbacks = []; // Store session feedbacks
        this.renderInitial();
        this.renderConsolidatedFeedback();
    },

    renderInitial() {
        this.container.innerHTML = `
            <div class="feedback-options">
                <button onclick="FeedbackSystem.openForm('agent_prompt')">Update Agent Prompt</button>
                <button onclick="FeedbackSystem.openForm('rubric')">Update Rubric</button>
                <button onclick="FeedbackSystem.openForm('question')">Clarify Question</button>
            </div>
            <div id="feedback-active-form"></div>
            <div id="feedback-history"></div>
        `;
        this.modal = document.getElementById('feedback-modal');
        this.modalTitle = document.getElementById('modal-title');
        this.modalContext = document.getElementById('modal-context');
        this.modalText = document.getElementById('modal-feedback-text');
        this.modalSubmit = document.getElementById('modal-submit-btn');
    },

    getFilePathForCategory(category) {
        if (category === 'rubric') return this.manifest.rubric?.file_path || '';
        if (category === 'agent_prompt') return this.manifest.prompts?.agent_system_source || '';
        if (category === 'question') return this.manifest.question?.file_path || '';
        return '';
    },

    openForm(category, target = '') {
        const titles = {
            'agent_prompt': 'Improve Agent Prompt',
            'rubric': 'Refine Rubric & Grading',
            'question': 'Clarify Task Question',
            'action': 'Specific Action Feedback'
        };

        this.modalTitle.textContent = titles[category];
        this.modalText.value = '';
        
        let contextHtml = `<div class="context-tag">${category.toUpperCase()}</div>`;
        
        const filePath = this.getFilePathForCategory(category);
        if (filePath) {
            contextHtml += `<div class="context-box"><strong>Target File:</strong> <code>${filePath}</code></div>`;
        }

        if (category === 'rubric' && target.startsWith('criterion:')) {
            const key = target.replace('criterion:', '');
            const latestEval = this.manifest.evaluations?.[this.manifest.evaluations.length - 1] || {};
            const catData = latestEval.categories?.[key] || {};
            const criterion = (this.manifest.rubric?.criteria || []).find(c => c.key === key) || {};
            
            contextHtml += `
                <div class="context-box">
                    <strong>Criterion:</strong> ${criterion.name}<br>
                    <strong>Current Score:</strong> ${catData.mean}/${criterion.max}<br>
                    <strong>Reasoning:</strong> <em>"${Array.isArray(catData.reasoning) ? catData.reasoning[0] : catData.reasoning}"</em>
                    ${catData.related_steps ? `<br><strong>Related Steps:</strong> #${catData.related_steps.join(', #')}` : ''}
                    ${catData.related_files ? `<br><strong>Related Files:</strong> ${catData.related_files.join(', ')}` : ''}
                </div>
            `;
        } else if (category === 'action' && target.startsWith('step:')) {
            const idx = parseInt(target.replace('step:', ''));
            const log = this.manifest.tool_log?.[idx] || {};
            contextHtml += `
                <div class="context-box">
                    <strong>Tool Call:</strong> ${log.name}<br>
                    <strong>Status:</strong> ${log.status}<br>
                    <strong>Preview:</strong> <code>${log.preview}</code>
                </div>
            `;
        } else {
            contextHtml += `<div class="context-box">Target: ${target || 'General'}</div>`;
        }

        this.modalContext.innerHTML = contextHtml;
        this.modalSubmit.onclick = () => this.submit(category, target);
        
        this.modal.classList.remove('hidden');
    },

    closeModal() {
        this.modal.classList.add('hidden');
    },

    submit(category, target) {
        const suggestion = this.modalText.value;

        if (!suggestion) {
            alert('Please provide a suggestion.');
            return;
        }

        const filePath = this.getFilePathForCategory(category);

        const feedback = {
            id: Date.now(),
            category,
            target,
            suggestion,
            timestamp: new Date().toISOString(),
            filePath: filePath || undefined
        };

        this.feedbacks.push(feedback);
        this.renderHistory();
        this.renderConsolidatedFeedback();
        this.closeModal();
        
        // Trigger optimization engine
        if (window.OptimizationEngine) {
            OptimizationEngine.onNewFeedback(feedback);
        }
    },

    submitDirect(category, target, suggestion) {
        if (!suggestion || suggestion.trim() === '') return;

        // Map 'action' feedback to agent_prompt file path if appropriate
        let filePath = this.getFilePathForCategory(category);
        if (category === 'action') {
            // Action feedback usually implies improving the agent prompt
            filePath = this.getFilePathForCategory('agent_prompt');
        }

        const feedback = {
            id: Date.now(),
            category,
            target,
            suggestion,
            timestamp: new Date().toISOString(),
            filePath: filePath || undefined
        };

        this.feedbacks.push(feedback);
        this.renderHistory();
        this.renderConsolidatedFeedback();
        
        // Trigger optimization engine
        if (window.OptimizationEngine) {
            OptimizationEngine.onNewFeedback(feedback);
        }
        
        console.log(`Direct feedback submitted for ${category} / ${target}`);
    },

    cancel() {
        this.closeModal();
    },

    renderConsolidatedFeedback() {
        const textArea = document.getElementById('consolidated-feedback-text');
        if (!textArea) return;

        if (this.feedbacks.length === 0) {
            textArea.value = '';
            return;
        }

        const consolidated = this.feedbacks.map(fb => {
            const date = new Date(fb.timestamp).toLocaleTimeString([], { hour: '2-digit', minute: '2-digit' });
            let header = `[${date}] ${fb.category.toUpperCase()}`;
            if (fb.target) header += ` (${fb.target})`;
            if (fb.filePath) header += `\nTarget File: ${fb.filePath}`;
            
            return `${header}\n${fb.suggestion}\n`;
        }).join('\n---\n\n');

        textArea.value = consolidated;
    },

    copyConsolidatedFeedback() {
        const textArea = document.getElementById('consolidated-feedback-text');
        if (!textArea || !textArea.value) return;

        textArea.select();
        document.execCommand('copy');
        
        // Visual feedback on the button
        const btn = document.querySelector('#consolidated-feedback-section .small-btn');
        const originalText = btn.textContent;
        btn.textContent = 'Copied!';
        btn.classList.add('success');
        setTimeout(() => {
            btn.textContent = originalText;
            btn.classList.remove('success');
        }, 2000);
    },

    renderHistory() {
        const historyDiv = document.getElementById('feedback-history');
        if (this.feedbacks.length === 0) {
            historyDiv.innerHTML = '';
            return;
        }

        historyDiv.innerHTML = '<h5>Recent Feedback</h5>';
        this.feedbacks.slice().reverse().forEach(fb => {
            const item = document.createElement('div');
            item.className = 'history-item';
            item.innerHTML = `
                <div class="badge">${fb.category}</div>
                <small>${fb.target}</small>
                <p>${fb.suggestion}</p>
            `;
            historyDiv.appendChild(item);
        });
    }
};

