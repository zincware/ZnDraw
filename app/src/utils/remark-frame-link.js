// remark-frame-link.js
import { visit, SKIP } from "unist-util-visit";

const frameRegex = /(^|\s)@(\d+)/g;

export function remarkFrameLink() {
  return (tree) => {
    visit(tree, "text", (node, index, parent) => {
      if (parent.type !== "paragraph") return;

      const text = node.value;
      const newNodes = [];
      let lastIndex = 0;
      let match;

      while ((match = frameRegex.exec(text)) !== null) {
        if (match.index > lastIndex) {
          newNodes.push({
            type: "text",
            value: text.substring(lastIndex, match.index),
          });
        }

        // preserve leading space
        if (match[1] === " ") {
          newNodes.push({ type: "text", value: " " });
        }

        const frameId = parseInt(match[2], 10);

        newNodes.push({
          type: "frameLink",
          data: {
            hName: "frameLink", // 👈 custom node instead of <a>
            hProperties: { frame: frameId },
          },
          children: [{ type: "text", value: `@${frameId}` }],
        });

        lastIndex = frameRegex.lastIndex;
      }

      if (newNodes.length === 0) return;

      if (lastIndex < text.length) {
        newNodes.push({ type: "text", value: text.substring(lastIndex) });
      }

      parent.children.splice(index, 1, ...newNodes);
      return [SKIP, index + newNodes.length];
    });
  };
}
