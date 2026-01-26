interface MatrixDisplayProps {
  data: number[][];
  label: string;
}

export default function MatrixDisplay({ data, label }: MatrixDisplayProps) {
  return (
    <div className="bg-gray-800/50 border border-gray-700/50 rounded-lg p-2 sm:p-3 backdrop-blur-sm">
      <div className="text-center font-semibold text-primary mb-2 text-sm sm:text-base">{label}</div>
      <div className="relative">
        {/* Left bracket */}
        <div className="absolute -left-2 sm:-left-2.5 top-0 bottom-0 w-1.5 sm:w-2 border-2 border-r-0 border-primary/50 rounded-l" />
        {/* Right bracket */}
        <div className="absolute -right-2 sm:-right-2.5 top-0 bottom-0 w-1.5 sm:w-2 border-2 border-l-0 border-primary/50 rounded-r" />

        <div className="grid grid-cols-3 gap-px bg-gray-700/50 border-2 border-primary/30 rounded">
          {data.flat().map((value, index) => (
            <div
              key={index}
              className="bg-gray-900/80 p-1.5 sm:p-2 text-center font-mono text-xs sm:text-sm min-w-[50px] sm:min-w-[70px] text-gray-200"
            >
              {value.toFixed(2)}
            </div>
          ))}
        </div>
      </div>
    </div>
  );
}
