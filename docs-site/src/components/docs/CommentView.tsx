import type { Comment, CommentDisplayPart } from '../../types/typedoc';

interface CommentViewProps {
  comment: Comment | undefined;
}

function renderParts(parts: CommentDisplayPart[] | undefined): string {
  if (!parts) return '';
  return parts.map(part => {
    if (part.kind === 'code') {
      return part.text;
    }
    return part.text;
  }).join('');
}

export default function CommentView({ comment }: CommentViewProps) {
  if (!comment) {
    return null;
  }

  const summary = renderParts(comment.summary);
  const examples = comment.blockTags?.filter(tag => tag.tag === '@example') || [];
  const returns = comment.blockTags?.find(tag => tag.tag === '@returns');

  return (
    <div className="comment text-gray-300">
      {summary && <p className="mb-4">{summary}</p>}

      {returns && (
        <div className="mt-4">
          <strong className="text-gray-200">Returns:</strong> {renderParts(returns.content)}
        </div>
      )}

      {examples.length > 0 && (
        <div className="mt-4">
          <strong className="text-gray-200">Example:</strong>
          {examples.map((example, i) => (
            <pre key={i} className="mt-2 bg-gray-800/50 border border-gray-700/50 p-4 rounded-lg overflow-x-auto">
              <code className="text-gray-200">{renderParts(example.content)}</code>
            </pre>
          ))}
        </div>
      )}
    </div>
  );
}
