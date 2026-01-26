import { ThreeBackground } from "./ThreeBackground";

export function Background() {
  return (
    <>
      {/* Gradient glow - top left corner */}
      <div
        className="pointer-events-none absolute -left-48 -top-48 z-0 h-[500px] w-[500px] rounded-full opacity-50 blur-3xl"
        style={{
          background: "radial-gradient(circle, #0d4a3a 0%, transparent 70%)",
        }}
      />

      {/* Gradient glow - bottom right corner */}
      <div
        className="pointer-events-none absolute -bottom-48 -right-48 z-0 h-[500px] w-[500px] rounded-full opacity-50 blur-3xl"
        style={{
          background: "radial-gradient(circle, #0d4a3a 0%, transparent 70%)",
        }}
      />

      {/* Gradient glow - center */}
      <div
        className="pointer-events-none absolute left-1/2 top-1/2 z-0 h-[600px] w-[600px] -translate-x-1/2 -translate-y-1/2 rounded-full opacity-30 blur-3xl"
        style={{
          background: "radial-gradient(circle, #0d4a3a 0%, transparent 70%)",
        }}
      />

      {/* Grid pattern with fade - bottom layer */}
      <div className="grid-pattern pointer-events-none absolute inset-0 z-0" />

      {/* Three.js 3D Scene */}
      <div className="absolute inset-0 z-[1]">
        <ThreeBackground />
      </div>

      {/* Scanlines effect - on top of 3D */}
      <div className="scanline pointer-events-none absolute inset-0 z-[2]" />

      {/* Noise texture - on top */}
      <div className="noise-overlay pointer-events-none absolute inset-0 z-[3]" />

      {/* Corner accents */}
      <div className="pointer-events-none absolute top-0 left-0 z-[4] h-32 w-32 border-t border-l border-primary/10" />
      <div className="pointer-events-none absolute top-0 right-0 z-[4] h-32 w-32 border-t border-r border-primary/10" />
      <div className="pointer-events-none absolute bottom-0 left-0 z-[4] h-32 w-32 border-b border-l border-primary/10" />
      <div className="pointer-events-none absolute bottom-0 right-0 z-[4] h-32 w-32 border-b border-r border-primary/10" />
    </>
  );
}
