import { useRef, useMemo, useEffect } from "react";
import { useFrame, useThree } from "@react-three/fiber";
import { OrbitControls, Text, Billboard, Html } from "@react-three/drei";
import * as THREE from "three";
import type { SimNode, SimEdge } from "./useForceSimulation";
import type { PackageId } from "./types";
import { PACKAGE_COLORS, ALL_PACKAGES } from "./types";

interface ForceGraphProps {
  nodes: SimNode[];
  edges: SimEdge[];
  onStep: () => void;
  isSettled: boolean;
  hoveredNode: SimNode | null;
  onHover: (node: SimNode | null) => void;
}

// Colors
const WASM_COLOR = new THREE.Color("#10b981");
const SHARED_COLOR = new THREE.Color("#f59e0b");
const BINDING_COLOR = new THREE.Color("#22c55e");

// Precompute package colors dynamically from ALL_PACKAGES
const packageColors: Record<PackageId, THREE.Color> = Object.fromEntries(
  ALL_PACKAGES.map((pkg) => [pkg, new THREE.Color(PACKAGE_COLORS[pkg])]),
) as Record<PackageId, THREE.Color>;

// Temp objects for instancing
const tempObject = new THREE.Object3D();
const tempColor = new THREE.Color();

function getNodeColor(node: SimNode): THREE.Color {
  if (node.isShared) return SHARED_COLOR;
  if (node.kind === "wasm") return WASM_COLOR;
  return packageColors[node.packageId] || new THREE.Color("#888888");
}

// Hover highlight color
const HOVER_COLOR = new THREE.Color("#ffffff");

// Search highlight color
const SEARCH_HIGHLIGHT_COLOR = new THREE.Color("#ffff00");

// Hover debounce delay in milliseconds
const HOVER_DEBOUNCE_MS = 150;

// Node instances component with hover support
function NodeInstances({
  nodes,
  hoveredNode,
  onHover,
}: {
  nodes: SimNode[];
  hoveredNode: SimNode | null;
  onHover: (node: SimNode | null) => void;
}) {
  const meshRef = useRef<THREE.InstancedMesh>(null);
  const hoveredIdRef = useRef<string | null>(null);
  const hoverTimeoutRef = useRef<number | null>(null);
  const { camera, raycaster, pointer } = useThree();

  // Create geometry for core sphere with enough segments for smooth shading
  const sphereGeometry = useMemo(() => new THREE.SphereGeometry(1, 32, 24), []);
  // Use MeshStandardMaterial for proper lighting - instanceColor is used automatically
  const sphereMaterial = useMemo(() => {
    const mat = new THREE.MeshStandardMaterial({
      roughness: 0.35,
      metalness: 0.15,
      toneMapped: false, // Required for bloom - allows colors > 1
    });
    return mat;
  }, []);

  // Initialize instance colors on first render
  useEffect(() => {
    if (!meshRef.current || nodes.length === 0) return;

    const mesh = meshRef.current;

    for (let i = 0; i < nodes.length; i++) {
      const node = nodes[i];
      tempColor.copy(getNodeColor(node));
      mesh.setColorAt(i, tempColor);
    }

    if (mesh.instanceColor) mesh.instanceColor.needsUpdate = true;
  }, [nodes]);

  // Sync React state to local ref
  useEffect(() => {
    hoveredIdRef.current = hoveredNode?.id ?? null;
  }, [hoveredNode]);

  // Update instances each frame and handle hover detection
  useFrame(() => {
    if (!meshRef.current || nodes.length === 0) return;

    const mesh = meshRef.current;

    // Update all instances first
    const activeHoveredId = hoveredIdRef.current;

    for (let i = 0; i < nodes.length; i++) {
      const node = nodes[i];
      const isHovered = activeHoveredId === node.id;
      const searchScore = node.searchScore || 0;
      const hasSearchQuery = nodes.some((n) => (n.searchScore || 0) > 0);

      // Calculate scale based on hover and search score
      let scale = node.radius;
      if (isHovered) {
        scale *= 1.5;
      } else if (hasSearchQuery && searchScore > 0) {
        scale *= 1 + searchScore * 1.5;
      } else if (hasSearchQuery) {
        scale *= 0.5;
      }

      // Position and scale for core sphere
      tempObject.position.set(node.x, node.y, node.z);
      tempObject.rotation.set(0, 0, 0);
      tempObject.scale.setScalar(scale);
      tempObject.updateMatrix();
      mesh.setMatrixAt(i, tempObject.matrix);

      // Color: hover > search highlight > normal
      // Multiply colors to make them emissive for bloom effect
      const emissiveBoost = 2.0; // Boost colors above 1.0 for bloom

      if (isHovered) {
        tempColor.copy(HOVER_COLOR).multiplyScalar(emissiveBoost * 1.5);
      } else if (hasSearchQuery && searchScore > 0) {
        const baseColor = getNodeColor(node);
        tempColor
          .copy(baseColor)
          .lerp(SEARCH_HIGHLIGHT_COLOR, searchScore * 0.6)
          .multiplyScalar(emissiveBoost);
      } else if (hasSearchQuery) {
        tempColor.copy(getNodeColor(node)).multiplyScalar(0.3);
      } else {
        tempColor.copy(getNodeColor(node)).multiplyScalar(emissiveBoost);
      }
      mesh.setColorAt(i, tempColor);
    }

    mesh.instanceMatrix.needsUpdate = true;
    if (mesh.instanceColor) mesh.instanceColor.needsUpdate = true;

    // Compute bounding sphere for raycasting
    mesh.computeBoundingSphere();

    // Raycast for hover detection (use core sphere)
    raycaster.setFromCamera(pointer, camera);
    const intersects = raycaster.intersectObject(mesh);

    let newHoveredId: string | null = null;
    let newHoveredNode: SimNode | null = null;

    if (intersects.length > 0 && intersects[0].instanceId !== undefined) {
      const instanceId = intersects[0].instanceId;
      if (instanceId < nodes.length) {
        newHoveredNode = nodes[instanceId];
        newHoveredId = newHoveredNode.id;
      }
    }

    const currentHoveredId = hoveredIdRef.current;

    if (newHoveredId !== currentHoveredId) {
      if (hoverTimeoutRef.current !== null) {
        clearTimeout(hoverTimeoutRef.current);
        hoverTimeoutRef.current = null;
      }

      if (newHoveredNode) {
        hoveredIdRef.current = newHoveredId;
        onHover(newHoveredNode);
      } else {
        hoverTimeoutRef.current = window.setTimeout(() => {
          hoveredIdRef.current = null;
          onHover(null);
          hoverTimeoutRef.current = null;
        }, HOVER_DEBOUNCE_MS);
      }
    }
  });

  if (nodes.length === 0) return null;

  return (
    <instancedMesh
      ref={meshRef}
      args={[sphereGeometry, sphereMaterial, nodes.length]}
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
      if (edge.edgeType === "binding") {
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
          const color =
            packageColors[source.packageId] || new THREE.Color("#555555");
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
    geo.setAttribute("position", new THREE.BufferAttribute(positions, 3));
    geo.setAttribute("color", new THREE.BufferAttribute(colors, 3));
    return geo;
  }, [containmentEdges.length]);

  const bindingGeometry = useMemo(() => {
    const geo = new THREE.BufferGeometry();
    const positions = new Float32Array(bindingEdges.length * 6);
    geo.setAttribute("position", new THREE.BufferAttribute(positions, 3));
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
            opacity={0.12}
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
            opacity={0.3}
            linewidth={1}
            dashSize={0.5}
            gapSize={0.3}
          />
        </lineSegments>
      )}
    </>
  );
}

// Individual label component that updates position
function NodeLabel({ node }: { node: SimNode }) {
  const groupRef = useRef<THREE.Group>(null);

  // Update position each frame
  useFrame(() => {
    if (groupRef.current) {
      groupRef.current.position.set(node.x, node.y + node.radius + 0.5, node.z);
    }
  });

  // Font size based on node type
  const fontSize =
    node.type === "package" ? 1.5 : node.type === "wasmModule" ? 0.8 : 0.4;

  // Get display label (first line only for packages)
  const displayLabel =
    node.type === "package"
      ? node.label.split("\n")[0]
      : node.label.replace(/^_+/, ""); // Remove leading underscores for functions

  return (
    <group
      ref={groupRef}
      position={[node.x, node.y + node.radius + 0.5, node.z]}
    >
      <Billboard follow lockX={false} lockY={false} lockZ={false}>
        <Text
          fontSize={fontSize}
          color="white"
          anchorX="center"
          anchorY="bottom"
          outlineWidth={0.05}
          outlineColor="#000000"
          maxWidth={node.type === "function" ? 10 : 20}
        >
          {displayLabel}
        </Text>
      </Billboard>
    </group>
  );
}

// Labels for all nodes
function NodeLabels({ nodes }: { nodes: SimNode[] }) {
  return (
    <>
      {nodes.map((node) => (
        <NodeLabel key={node.id} node={node} />
      ))}
    </>
  );
}

// Center sun component - glowing sphere that acts as light source
function CenterSun() {
  const sunRef = useRef<THREE.Mesh>(null);

  // Subtle rotation animation
  useFrame((_, delta) => {
    if (sunRef.current) {
      sunRef.current.rotation.y += delta * 0.1;
    }
  });

  return (
    <group position={[0, 0, 0]}>
      {/* Point light emanating from the sun - high intensity, no distance limit */}
      <pointLight intensity={400} decay={1.1} color="#fff8e8" />

      {/* Glowing sun sphere */}
      <mesh ref={sunRef}>
        <sphereGeometry args={[8, 32, 32]} />
        <meshBasicMaterial color="#ffcc66" toneMapped={false} />
      </mesh>
    </group>
  );
}

// Background stars - fixed points in space
function BackgroundStars() {
  const pointsRef = useRef<THREE.Points>(null);

  const { geometry, material } = useMemo(() => {
    const starCount = 2000;
    const positions = new Float32Array(starCount * 3);
    const sizes = new Float32Array(starCount);

    // Distribute stars on a large sphere around the scene
    const radius = 500;
    for (let i = 0; i < starCount; i++) {
      // Random spherical distribution
      const theta = Math.random() * Math.PI * 2;
      const phi = Math.acos(2 * Math.random() - 1);

      positions[i * 3] = radius * Math.sin(phi) * Math.cos(theta);
      positions[i * 3 + 1] = radius * Math.sin(phi) * Math.sin(theta);
      positions[i * 3 + 2] = radius * Math.cos(phi);

      // Vary star sizes slightly
      sizes[i] = 0.5 + Math.random() * 1.5;
    }

    const geo = new THREE.BufferGeometry();
    geo.setAttribute("position", new THREE.BufferAttribute(positions, 3));
    geo.setAttribute("size", new THREE.BufferAttribute(sizes, 1));

    const mat = new THREE.PointsMaterial({
      color: 0xffffff,
      size: 1,
      sizeAttenuation: true,
      transparent: true,
      opacity: 0.8,
    });

    return { geometry: geo, material: mat };
  }, []);

  return <points ref={pointsRef} geometry={geometry} material={material} />;
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

export function ForceGraph({
  nodes,
  edges,
  onStep,
  isSettled,
  hoveredNode,
  onHover,
}: ForceGraphProps) {
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

      {/* Background stars */}
      <BackgroundStars />

      {/* Center Sun - light source */}
      <CenterSun />

      {/* Minimal ambient so nodes aren't completely dark on back side */}
      <ambientLight intensity={0.1} />

      {/* Nodes */}
      <NodeInstances
        nodes={nodes}
        hoveredNode={hoveredNode}
        onHover={onHover}
      />

      {/* Edges */}
      <EdgeLines nodes={nodes} edges={edges} />

      {/* Labels */}
      <NodeLabels nodes={nodes} />

      {/* Tooltip for hovered node */}
      {hoveredNode && <HoverTooltip node={hoveredNode} />}
    </>
  );
}

// HoverTooltip - positioned based on screen location to avoid covering the node
function HoverTooltip({ node }: { node: SimNode }) {
  const { camera } = useThree();
  const groupRef = useRef<THREE.Group>(null);

  // Helper to calculate offset scaled by camera distance
  const calculateOffset = () => {
    const nodePos = new THREE.Vector3(node.x, node.y, node.z);

    // Calculate distance from camera to node
    const cameraDistance = camera.position.distanceTo(nodePos);

    // Scale offset based on distance (larger offset when zoomed out)
    // Base offset of 4 at distance 80, scales proportionally
    const scaleFactor = cameraDistance / 80;
    const baseOffsetX = 6 * scaleFactor;
    const baseOffsetY = 5 * scaleFactor;

    // Project to get screen position
    nodePos.project(camera);
    const screenX = (nodePos.x + 1) / 2;
    const screenY = (nodePos.y + 1) / 2;

    // Determine direction based on screen position
    const offsetX = screenX > 0.5 ? -baseOffsetX : baseOffsetX;
    const offsetY = screenY > 0.5 ? -baseOffsetY : baseOffsetY;

    return { offsetX, offsetY, screenX, screenY };
  };

  // Calculate screen position and determine tooltip placement
  const tooltipPosition = useMemo(() => {
    const { offsetX, offsetY } = calculateOffset();
    return {
      x: node.x + offsetX,
      y: node.y + offsetY,
      z: node.z,
    };
  }, [node.x, node.y, node.z, camera]);

  // Update position each frame to follow node and adjust for camera movement
  useFrame(() => {
    if (groupRef.current) {
      const { offsetX, offsetY } = calculateOffset();
      groupRef.current.position.set(node.x + offsetX, node.y + offsetY, node.z);
    }
  });

  const getTypeLabel = () => {
    if (node.type === "package") return "Package";
    if (node.type === "wasmModule") return "WASM Module";
    if (node.kind === "wasm") return "WASM Function";
    if (node.kind === "typescript") return "TypeScript Export";
    return "Function";
  };

  const getDescription = () => {
    // Use actual description if available
    if (node.description) {
      return node.description;
    }

    // Fallback descriptions for nodes without docs
    if (node.type === "package") {
      return `Root package containing WASM modules and TypeScript exports`;
    }
    if (node.type === "wasmModule") {
      return `WebAssembly module compiled from native code`;
    }
    if (node.kind === "wasm") {
      const name = node.label.replace(/^_+/, "");
      if (name.includes("create")) return "Creates a new instance";
      if (name.includes("free")) return "Frees allocated memory";
      if (name.includes("get")) return "Retrieves a value or property";
      if (name.includes("set")) return "Sets a value or property";
      if (name.includes("add")) return "Adds or combines values";
      if (name.includes("multiply") || name.includes("mul"))
        return "Multiplies values";
      if (name.includes("divide") || name.includes("div"))
        return "Divides values";
      if (name.includes("solve")) return "Solves a mathematical problem";
      if (name.includes("factor")) return "Performs factorization";
      if (name.includes("eigen")) return "Computes eigenvalues/eigenvectors";
      if (name.includes("sort")) return "Sorts elements";
      if (name.includes("search")) return "Searches for elements";
      return "Low-level WASM function";
    }
    if (node.kind === "typescript") {
      return "High-level TypeScript API function";
    }
    return "";
  };

  return (
    <group
      ref={groupRef}
      position={[tooltipPosition.x, tooltipPosition.y, tooltipPosition.z]}
    >
      <Html center style={{ pointerEvents: "none" }}>
        <div
          style={{
            background: "rgba(15, 15, 26, 0.95)",
            border: "1px solid rgba(255, 255, 255, 0.2)",
            borderRadius: "8px",
            padding: "12px 16px",
            minWidth: "200px",
            maxWidth: "300px",
            boxShadow: "0 4px 20px rgba(0, 0, 0, 0.5)",
          }}
        >
          <div
            style={{
              fontSize: "10px",
              textTransform: "uppercase",
              letterSpacing: "0.5px",
              color: node.color,
              marginBottom: "4px",
            }}
          >
            {getTypeLabel()}
          </div>
          <div
            style={{
              fontSize: "14px",
              fontWeight: "bold",
              color: "white",
              marginBottom: "8px",
              fontFamily: "monospace",
              wordBreak: "break-all",
            }}
          >
            {node.label.replace(/^_+/, "")}
          </div>
          <div
            style={{
              fontSize: "12px",
              color: "rgba(255, 255, 255, 0.7)",
              lineHeight: "1.4",
            }}
          >
            {getDescription()}
          </div>
          {node.isShared && (
            <div
              style={{
                marginTop: "8px",
                fontSize: "11px",
                color: "#f59e0b",
                display: "flex",
                alignItems: "center",
                gap: "4px",
              }}
            >
              ⚠ Shared across multiple packages
            </div>
          )}
          <div
            style={{
              marginTop: "8px",
              fontSize: "11px",
              color: "rgba(255, 255, 255, 0.5)",
            }}
          >
            Package: {node.packageId}
          </div>
        </div>
      </Html>
    </group>
  );
}
