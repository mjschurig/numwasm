import { Preload } from "@react-three/drei";
import { Canvas } from "@react-three/fiber";
import { EnergySphere } from "./EnergySphere";
import { Particles } from "./Particles";

function Scene() {
  return (
    <>
      {/* Floating ambient particles */}
      <Particles count={150} />

      {/* Energy sphere - positioned to top right */}
      <EnergySphere position={[5, 4, -4]} scale={5} />

      <Preload all />
    </>
  );
}

export function ThreeBackground() {
  return (
    <Canvas
      camera={{ position: [0, 0, 5], fov: 75 }}
      dpr={[1, 1.5]} // Limit pixel ratio for performance
      gl={{
        antialias: true,
        alpha: true, // Transparent background
        powerPreference: "high-performance",
        preserveDrawingBuffer: true,
      }}
      style={{
        position: "absolute",
        top: 0,
        left: 0,
        width: "100%",
        height: "100%",
        pointerEvents: "none",
        zIndex: 0,
      }}
    >
      <Scene />
    </Canvas>
  );
}
