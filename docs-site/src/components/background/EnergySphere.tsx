import { useFrame } from "@react-three/fiber";
import { useMemo, useRef } from "react";
import * as THREE from "three";

// Very dark, subtle cyan for background
const CYAN_DIM = "#0d4a3a";

interface EnergySphereProps {
  position?: [number, number, number];
  scale?: number;
}

export function EnergySphere({
  position = [0, 0, 0],
  scale = 1,
}: EnergySphereProps) {
  const groupRef = useRef<THREE.Group>(null);
  const sphereRef = useRef<THREE.LineSegments>(null);

  // Create geometries
  const sphereGeo = useMemo(() => {
    // Wireframe sphere - icosahedron for geometric look
    const ico = new THREE.IcosahedronGeometry(1.5, 1);
    return new THREE.EdgesGeometry(ico, 1);
  }, []);

  useFrame((state) => {
    const time = state.clock.elapsedTime;

    if (groupRef.current) {
      // Slow rotation
      groupRef.current.rotation.y = time * 0.06;
      groupRef.current.rotation.x = time * 0.02;
    }

    // Sphere gentle pulsing
    if (sphereRef.current) {
      const pulse = 1 + Math.sin(time * 1.5) * 0.02;
      sphereRef.current.scale.setScalar(pulse);
    }
  });

  return (
    <group ref={groupRef} position={position} scale={scale}>
      {/* Core wireframe sphere */}
      <lineSegments ref={sphereRef} geometry={sphereGeo}>
        <lineBasicMaterial color={CYAN_DIM} transparent opacity={0.35} />
      </lineSegments>
    </group>
  );
}
