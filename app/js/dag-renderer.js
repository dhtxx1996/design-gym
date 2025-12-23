const DagRenderer = {
    init(manifest) {
        this.canvas = document.getElementById('dag-canvas');
        this.dag = manifest.dag || { nodes: [], edges: [] };
        this.render();
    },

    render() {
        if (!this.dag.nodes.length) {
            this.canvas.innerHTML = '<div class="empty">No DAG data available.</div>';
            return;
        }

        const width = this.canvas.clientWidth || 800;
        const height = this.canvas.clientHeight || 600;

        // Simple Layout Algorithm (Vertical Hierarchical)
        const layout = this.calculateLayout(width, height);
        
        let svgHtml = `<svg width="${width}" height="${height}" viewBox="0 0 ${width} ${height}">`;
        
        // Render edges first (so they are behind nodes)
        this.dag.edges.forEach(edge => {
            const start = layout[edge.source];
            const end = layout[edge.target];
            if (start && end) {
                svgHtml += `
                    <line x1="${start.x}" y1="${start.y}" x2="${end.x}" y2="${end.y}" 
                          stroke="#adb5bd" stroke-width="2" marker-end="url(#arrowhead)" />
                `;
            }
        });

        // Render nodes
        this.dag.nodes.forEach(node => {
            const pos = layout[node.id];
            if (!pos) return;

            const isFunction = node.type === 'function';
            const color = isFunction ? '#e7f5ff' : '#fff9db';
            const stroke = isFunction ? '#228be6' : '#fab005';
            
            // Draw rectangle instead of circle
            const padding = 8;
            const fontSize = 10;
            const textWidth = node.label.length * 6; // Rough estimation
            const rectWidth = Math.max(textWidth + padding * 2, 60);
            const rectHeight = 24;

            svgHtml += `
                <g class="dag-node" data-id="${node.id}" transform="translate(${pos.x - rectWidth/2}, ${pos.y - rectHeight/2})"
                   onclick="DagRenderer.onNodeClick('${node.id}')">
                    <rect width="${rectWidth}" height="${rectHeight}" rx="6" ry="6" 
                          fill="${color}" stroke="${stroke}" stroke-width="2" />
                    <text x="${rectWidth/2}" y="${rectHeight/2 + 4}" text-anchor="middle" font-size="${fontSize}" fill="#495057">${node.label}</text>
                </g>
            `;
        });

        // Arrowhead definition
        svgHtml += `
            <defs>
                <marker id="arrowhead" markerWidth="10" markerHeight="7" refX="10" refY="3.5" orient="auto">
                    <polygon points="0 0, 10 3.5, 0 7" fill="#adb5bd" />
                </marker>
            </defs>
        </svg>`;

        this.canvas.innerHTML = svgHtml;
    },

    calculateLayout(width, height) {
        const layout = {};
        const nodes = this.dag.nodes;
        const edges = this.dag.edges;

        // Group nodes by "layers" based on incoming edges
        const inDegrees = {};
        nodes.forEach(n => inDegrees[n.id] = 0);
        edges.forEach(e => inDegrees[e.target]++);

        let currentLayer = nodes.filter(n => inDegrees[n.id] === 0).map(n => n.id);
        let layerIdx = 0;
        const layers = [];

        while (currentLayer.length > 0) {
            layers.push(currentLayer);
            const nextLayer = [];
            currentLayer.forEach(nodeId => {
                edges.filter(e => e.source === nodeId).forEach(e => {
                    inDegrees[e.target]--;
                    if (inDegrees[e.target] === 0) {
                        nextLayer.push(e.target);
                    }
                });
            });
            currentLayer = nextLayer;
            layerIdx++;
        }

        // Assign coordinates based on layers
        const layerHeight = height / (layers.length + 1);
        layers.forEach((layer, lIdx) => {
            const y = (lIdx + 1) * layerHeight;
            const layerWidth = width / (layer.length + 1);
            layer.forEach((nodeId, nIdx) => {
                layout[nodeId] = {
                    x: (nIdx + 1) * layerWidth,
                    y: y
                };
            });
        });

        return layout;
    },

    onNodeClick(nodeId) {
        console.log('Node clicked:', nodeId);
        
        // Clear previous highlights
        this.clearHighlights();
        this.highlightNode(nodeId);

        // If it's a tool call (function), highlight it in the trace viewer
        // Extract tool call index from node id (e.g., run_python_1 -> idx 1)
        const match = nodeId.match(/_(\d+)$/);
        if (match) {
            const idx = parseInt(match[1]);
            console.log('Detected step index from node:', idx);
            if (window.TraceViewer) {
                window.TraceViewer.highlight(idx);
            }
            if (window.RubricViewer) {
                window.RubricViewer.highlightForStep(idx);
            }
        } else if (nodeId.startsWith('file:')) {
            const fileName = nodeId.replace('file:', '');
            if (window.RubricViewer) {
                window.RubricViewer.highlightForFile(fileName);
            }
        }
    },

    highlightForStep(idx) {
        const nodeIdPattern = `_${idx}`;
        const nodes = this.canvas.querySelectorAll('.dag-node');
        nodes.forEach(node => {
            node.classList.remove('highlighted');
            if (node.dataset.id.endsWith(nodeIdPattern)) {
                node.classList.add('highlighted');
            }
        });
    },

    highlightNode(nodeId) {
        const nodes = this.canvas.querySelectorAll('.dag-node');
        nodes.forEach(node => {
            if (node.dataset.id === nodeId) {
                node.classList.add('highlighted');
            }
        });
    },

    clearHighlights() {
        const nodes = this.canvas.querySelectorAll('.dag-node');
        nodes.forEach(node => node.classList.remove('highlighted'));
    }
};

