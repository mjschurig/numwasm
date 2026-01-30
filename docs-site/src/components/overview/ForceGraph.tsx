import { useRef, useMemo, useEffect } from 'react';
import { useFrame, useThree } from '@react-three/fiber';
import { OrbitControls, Html } from '@react-three/drei';
import * as THREE from 'three';
import type { SimNode, SimEdge } from './useForceSimulation';
import type { PackageId } from './types';
import { PACKAGE_COLORS } from './types';

interface ForceGraphProps {
  nodes: SimNode[];
  edges: SimEdge[];
  onStep: () => void;
  isSettled: boolean;
}

// Colors
const WASM_COLOR = new THREE.Color('#10b981');
const SHARED_COLOR = new THREE.Color('#f59e0b');
const BINDING_COLOR = new THREE.Color('#22c55e');

// Precompute package colors
const packageColors: Record<PackageId, THREE.Color> = {
  numwasm: new THREE.Color(PACKAGE_COLORS.numwasm),
  sciwasm: new THREE.Color(PACKAGE_COLORS.sciwasm),
  symwasm: new THREE.Color(PACKAGE_COLORS.symwasm),
};

// Temp objects for instancing
const tempObject = new THREE.Object3D();
const tempColor = new THREE.Color();

function getNodeColor(node: SimNode): THREE.Color {
  if (node.isShared) return SHARED_COLOR;
  if (node.kind === 'wasm') return WASM_COLOR;
  return packageColors[node.packageId] || new THREE.Color('#888888');
}

// Node instances component
function NodeInstances({ nodes }: { nodes: SimNode[] }) {
  const meshRef = useRef<THREE.InstancedMesh>(null);

  // Create geometry and material
  const geometry = useMemo(() => new THREE.SphereGeometry(1, 16, 12), []);
  const material = useMemo(
    () =>
      new THREE.MeshBasicMaterial({
        vertexColors: true,
        transparent: true,
        opacity: 0.9,
      }),
    []
  );

  // Update instances each frame
  useFrame(() => {
    if (!meshRef.current || nodes.length === 0) return;

    const mesh = meshRef.current;

    for (let i = 0; i < nodes.length; i++) {
      const node = nodes[i];

      // Position and scale
      tempObject.position.set(node.x, node.y, node.z);
      tempObject.scale.setScalar(node.radius);
      tempObject.updateMatrix();
      mesh.setMatrixAt(i, tempObject.matrix);

      // Color
      tempColor.copy(getNodeColor(node));
      mesh.setColorAt(i, tempColor);
    }

    mesh.instanceMatrix.needsUpdate = true;
    if (mesh.instanceColor) mesh.instanceColor.needsUpdate = true;
  });

  if (nodes.length === 0) return null;

  return (
    <instancedMesh
      ref={meshRef}
      args={[geometry, material, nodes.length]}
      frustumCulled={false}
    />
  );
}

// Edge lines component
function EdgeLines({ nodes, edges }: { nodes: SimNode[]; edges: SimEdge[] }) {
  const containmentRef = useRef<THREE.LineSegments>(null);
  const bindingRef = useRef<THREE.LineSegments>(null);

  // Create node lookup map
  const nodeMap = useMemo(() => {
    const map = new Map<string, SimNode>();
    for (const node of nodes) {
      map.set(node.id, node);
    }
    return map;
  }, [nodes]);

  // Separate edges by type
  const { containmentEdges, bindingEdges } = useMemo(() => {
    const containment: SimEdge[] = [];
    const binding: SimEdge[] = [];
    for (const edge of edges) {
      if (edge.edgeType === 'binding') {
        binding.push(edge);
      } else {
        containment.push(edge);
      }
    }
    return { containmentEdges: containment, bindingEdges: binding };
  }, [edges]);

  // Update edge positions each frame
  useFrame(() => {
    // Update containment edges
    if (containmentRef.current && containmentEdges.length > 0) {
      const positions = containmentRef.current.geometry.attributes.position
        .array as Float32Array;
      const colors = containmentRef.current.geometry.attributes.color
        .array as Float32Array;

      for (let i = 0; i < containmentEdges.length; i++) {
        const edge = containmentEdges[i];
        const source = nodeMap.get(edge.source);
        const target = nodeMap.get(edge.target);

        if (source && target) {
          const idx = i * 6;
          positions[idx] = source.x;
          positions[idx + 1] = source.y;
          positions[idx + 2] = source.z;
          positions[idx + 3] = target.x;
          positions[idx + 4] = target.y;
          positions[idx + 5] = target.z;

          // Color from source node's package
          const color = packageColors[source.packageId] || new THREE.Color('#555555');
          colors[idx] = color.r;
          colors[idx + 1] = color.g;
          colors[idx + 2] = color.b;
          colors[idx + 3] = color.r;
          colors[idx + 4] = color.g;
          colors[idx + 5] = color.b;
        }
      }

      containmentRef.current.geometry.attributes.position.needsUpdate = true;
      containmentRef.current.geometry.attributes.color.needsUpdate = true;
    }

    // Update binding edges
    if (bindingRef.current && bindingEdges.length > 0) {
      const positions = bindingRef.current.geometry.attributes.position
        .array as Float32Array;

      for (let i = 0; i < bindingEdges.length; i++) {
        const edge = bindingEdges[i];
        const source = nodeMap.get(edge.source);
        const target = nodeMap.get(edge.target);

        if (source && target) {
          const idx = i * 6;
          positions[idx] = source.x;
          positions[idx + 1] = source.y;
          positions[idx + 2] = source.z;
          positions[idx + 3] = target.x;
          positions[idx + 4] = target.y;
          positions[idx + 5] = target.z;
        }
      }

      bindingRef.current.geometry.attributes.position.needsUpdate = true;
    }
  });

  // Create geometries
  const containmentGeometry = useMemo(() => {
    const geo = new THREE.BufferGeometry();
    const positions = new Float32Array(containmentEdges.length * 6);
    const colors = new Float32Array(containmentEdges.length * 6);
    geo.setAttribute('position', new THREE.BufferAttribute(positions, 3));
    geo.setAttribute('color', new THREE.BufferAttribute(colors, 3));
    return geo;
  }, [containmentEdges.length]);

  const bindingGeometry = useMemo(() => {
    const geo = new THREE.BufferGeometry();
    const positions = new Float32Array(bindingEdges.length * 6);
    geo.setAttribute('position', new THREE.BufferAttribute(positions, 3));
    return geo;
  }, [bindingEdges.length]);

  return (
    <>
      {/* Containment edges - solid, colored by package */}
      {containmentEdges.length > 0 && (
        <lineSegments ref={containmentRef} geometry={containmentGeometry}>
          <lineBasicMaterial
            vertexColors
            transparent
            opacity={0.3}
            linewidth={1}
          />
        </lineSegments>
      )}

      {/* Binding edges - dashed green */}
      {bindingEdges.length > 0 && (
        <lineSegments ref={bindingRef} geometry={bindingGeometry}>
          <lineDashedMaterial
            color={BINDING_COLOR}
            transparent
            opacity={0.7}
            linewidth={2}
            dashSize={0.5}
            gapSize={0.3}
          />
        </lineSegments>
      )}
    </>
  );
}

// Labels for package nodes only (to avoid clutter)
function NodeLabels({ nodes }: { nodes: SimNode[] }) {
  const packageNodes = useMemo(
    () => nodes.filter((n) => n.type === 'package'),
    [nodes]
  );

  return (
    <>
      {packageNodes.map((node) => (
        <Html
          key={node.id}
          position={[node.x, node.y + node.radius + 1, node.z]}
          center
          style={{
            color: 'white',
            fontSize: '14px',
            fontWeight: 'bold',
            textShadow: '0 0 4px rgba(0,0,0,0.8)',
            whiteSpace: 'nowrap',
            pointerEvents: 'none',
          }}
        >
          {node.label.split('\n')[0]}
        </Html>
      ))}
    </>
  );
}

// Camera setup
function CameraSetup() {
  const { camera } = useThree();

  useEffect(() => {
    camera.position.set(0, 0, 80);
    camera.lookAt(0, 0, 0);
  }, [camera]);

  return null;
}

export function ForceGraph({ nodes, edges, onStep, isSettled }: ForceGraphProps) {
  // Run simulation step each frame until settled
  useFrame(() => {
    if (!isSettled) {
      onStep();
    }
  });

  return (
    <>
      <CameraSetup />
      <OrbitControls
        enableDamping
        dampingFactor={0.1}
        rotateSpeed={0.5}
        zoomSpeed={1.2}
        panSpeed={0.8}
        minDistance={10}
        maxDistance={200}
      />

      {/* Ambient light for visibility */}
      <ambientLight intensity={1} />

      {/* Nodes */}
      <NodeInstances nodes={nodes} />

      {/* Edges */}
      <EdgeLines nodes={nodes} edges={edges} />

      {/* Labels */}
      <NodeLabels nodes={nodes} />
    </>
  );
}
