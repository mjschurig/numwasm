import { useFrame } from "@react-three/fiber";
import { useEffect, useMemo, useRef } from "react";
import * as THREE from "three";

// Brand color - cyan
const CYAN = "#2dd4a8";

// Create a circular texture for rounded particles
function createCircleTexture(): THREE.Texture {
  const canvas = document.createElement("canvas");
  canvas.width = 64;
  canvas.height = 64;
  const ctx = canvas.getContext("2d");
  if (ctx) {
    const gradient = ctx.createRadialGradient(32, 32, 0, 32, 32, 32);
    gradient.addColorStop(0, "rgba(255, 255, 255, 1)");
    gradient.addColorStop(0.5, "rgba(255, 255, 255, 0.8)");
    gradient.addColorStop(1, "rgba(255, 255, 255, 0)");
    ctx.fillStyle = gradient;
    ctx.fillRect(0, 0, 64, 64);
  }
  const texture = new THREE.CanvasTexture(canvas);
  return texture;
}

interface ParticlesProps {
  count?: number;
}

// Floating ambient particles with varied sizes
export function Particles({ count = 80 }: ParticlesProps) {
  const smallRef = useRef<THREE.Points>(null);
  const mediumRef = useRef<THREE.Points>(null);
  const largeRef = useRef<THREE.Points>(null);
  const smallMatRef = useRef<THREE.PointsMaterial>(null);
  const mediumMatRef = useRef<THREE.PointsMaterial>(null);
  const largeMatRef = useRef<THREE.PointsMaterial>(null);

  // Create circle texture on mount
  useEffect(() => {
    const texture = createCircleTexture();
    if (smallMatRef.current) {
      smallMatRef.current.map = texture;
      smallMatRef.current.needsUpdate = true;
    }
    if (mediumMatRef.current) {
      mediumMatRef.current.map = texture.clone();
      mediumMatRef.current.needsUpdate = true;
    }
    if (largeMatRef.current) {
      largeMatRef.current.map = texture.clone();
      largeMatRef.current.needsUpdate = true;
    }
  }, []);

  // Generate particle data for different sizes
  const { small, medium, large } = useMemo(() => {
    const smallCount = Math.floor(count * 0.7);
    const mediumCount = Math.floor(count * 0.2);
    const largeCount = Math.floor(count * 0.1);

    const cyanColor = new THREE.Color(CYAN);

    const createParticles = (num: number) => {
      const positions = new Float32Array(num * 3);
      const colors = new Float32Array(num * 3);
      for (let i = 0; i < num; i++) {
        positions[i * 3] = (Math.random() - 0.5) * 25;
        positions[i * 3 + 1] = (Math.random() - 0.5) * 15;
        positions[i * 3 + 2] = -8 + (Math.random() - 0.5) * 8;
        colors[i * 3] = cyanColor.r;
        colors[i * 3 + 1] = cyanColor.g;
        colors[i * 3 + 2] = cyanColor.b;
      }
      return { positions, colors, count: num };
    };

    return {
      small: createParticles(smallCount),
      medium: createParticles(mediumCount),
      large: createParticles(largeCount),
    };
  }, [count]);

  // Floating animation with drift
  useFrame((state) => {
    const time = state.clock.elapsedTime;

    const animateParticles = (
      ref: React.RefObject<THREE.Points | null>,
      particleCount: number,
      speedMult: number,
    ) => {
      if (ref.current) {
        const positions = ref.current.geometry.attributes.position
          .array as Float32Array;
        for (let i = 0; i < particleCount; i++) {
          const offset = i * 0.15;
          const speed = (0.5 + (i % 5) * 0.1) * speedMult;
          positions[i * 3] += Math.sin(time * speed + offset) * 0.003;
          positions[i * 3 + 1] += Math.cos(time * speed * 0.8 + offset) * 0.004;
          positions[i * 3 + 2] +=
            Math.sin(time * speed * 0.5 + offset * 2) * 0.001;
        }
        ref.current.geometry.attributes.position.needsUpdate = true;
      }
    };

    animateParticles(smallRef, small.count, 1);
    animateParticles(mediumRef, medium.count, 0.7);
    animateParticles(largeRef, large.count, 0.4);
  });

  return (
    <>
      {/* Small particles */}
      <points ref={smallRef}>
        <bufferGeometry>
          <bufferAttribute
            attach="attributes-position"
            args={[small.positions, 3]}
          />
          <bufferAttribute attach="attributes-color" args={[small.colors, 3]} />
        </bufferGeometry>
        <pointsMaterial
          ref={smallMatRef}
          size={0.03}
          vertexColors
          transparent
          opacity={0.4}
          sizeAttenuation
          depthWrite={false}
          blending={THREE.AdditiveBlending}
        />
      </points>

      {/* Medium particles */}
      <points ref={mediumRef}>
        <bufferGeometry>
          <bufferAttribute
            attach="attributes-position"
            args={[medium.positions, 3]}
          />
          <bufferAttribute
            attach="attributes-color"
            args={[medium.colors, 3]}
          />
        </bufferGeometry>
        <pointsMaterial
          ref={mediumMatRef}
          size={0.08}
          vertexColors
          transparent
          opacity={0.5}
          sizeAttenuation
          depthWrite={false}
          blending={THREE.AdditiveBlending}
        />
      </points>

      {/* Large particles */}
      <points ref={largeRef}>
        <bufferGeometry>
          <bufferAttribute
            attach="attributes-position"
            args={[large.positions, 3]}
          />
          <bufferAttribute attach="attributes-color" args={[large.colors, 3]} />
        </bufferGeometry>
        <pointsMaterial
          ref={largeMatRef}
          size={0.15}
          vertexColors
          transparent
          opacity={0.6}
          sizeAttenuation
          depthWrite={false}
          blending={THREE.AdditiveBlending}
        />
      </points>
    </>
  );
}
