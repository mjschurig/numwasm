import { useRef, useCallback, useMemo } from "react";
import type { CyElement, CyNodeData, CyEdgeData, PackageId } from "./types";

export interface SimNode {
  id: string;
  x: number;
  y: number;
  z: number;
  vx: number;
  vy: number;
  vz: number;
  type: "package" | "wasmModule" | "function";
  packageId: PackageId;
  color: string;
  label: string;
  mass: number;
  radius: number;
  kind?: "wasm" | "typescript";
  isShared?: boolean;
  description?: string;
  searchScore?: number; // 0-1, how well this node matches the search query
}

export interface SimEdge {
  source: string;
  target: string;
  edgeType: "containment" | "binding";
  color: string;
}

interface ForceSimulationResult {
  nodes: SimNode[];
  edges: SimEdge[];
  step: () => void;
  isSettled: boolean;
}

// Force parameters
const FORCES = {
  repulsion: 500,
  attraction: 0.1,
  bindingAttraction: 2, // Stronger attraction for TS→WASM binding edges
  centerPull: 0.001,
  clusterPull: 50,
  packageDistance: 120, // Fixed distance from center for all packages
  damping: 0.95, // Higher damping for faster settling
  maxVelocity: 5,
  minDistance: 0.5,
  // Simulated annealing - reduce forces over time
  initialAlpha: 0.5,
  alphaDecay: 0.01, // How fast alpha decreases per iteration
  alphaMin: 0.005, // Stop when alpha reaches this
  velocityDecay: 0.2, // Additional velocity decay per frame
};

// Node sizes by type
const NODE_RADIUS: Record<string, number> = {
  package: 2.0,
  wasmModule: 1.2,
  function: 0.5,
};

const NODE_MASS: Record<string, number> = {
  package: 20,
  wasmModule: 5,
  function: 1,
};

function isNodeElement(el: CyElement): el is CyElement & { data: CyNodeData } {
  return !("source" in el.data);
}

function isEdgeElement(el: CyElement): el is CyElement & { data: CyEdgeData } {
  return "source" in el.data;
}

/**
 * Calculate fuzzy string similarity between query and target
 * Returns a score between 0 and 1
 */
function calculateSearchScore(query: string, target: string): number {
  if (!query) return 0;

  const q = query.toLowerCase();
  const t = target.toLowerCase().replace(/^_+/, ""); // Remove leading underscores

  // Exact match
  if (t === q) return 1;

  // Starts with query (high score)
  if (t.startsWith(q)) return 0.9;

  // Contains query as substring
  if (t.includes(q)) {
    // Score based on how early the match appears and relative length
    const index = t.indexOf(q);
    const positionScore = 1 - (index / t.length) * 0.3;
    const lengthScore = q.length / t.length;
    return Math.min(0.85, 0.5 + positionScore * 0.2 + lengthScore * 0.15);
  }

  // Fuzzy matching: check if all characters appear in order
  let queryIdx = 0;
  let matchCount = 0;
  let consecutiveBonus = 0;
  let lastMatchIdx = -2;

  for (let i = 0; i < t.length && queryIdx < q.length; i++) {
    if (t[i] === q[queryIdx]) {
      matchCount++;
      if (i === lastMatchIdx + 1) {
        consecutiveBonus += 0.1;
      }
      lastMatchIdx = i;
      queryIdx++;
    }
  }

  if (queryIdx === q.length) {
    // All query characters found in order
    const baseScore = matchCount / t.length;
    return Math.min(0.6, baseScore * 0.4 + consecutiveBonus);
  }

  // Partial character match (some characters found)
  if (matchCount > 0) {
    return (matchCount / q.length) * 0.2;
  }

  return 0;
}

export function useForceSimulation(
  elements: CyElement[],
  searchQuery: string = "",
): ForceSimulationResult {
  const nodesRef = useRef<SimNode[]>([]);
  const edgesRef = useRef<SimEdge[]>([]);
  const nodeMapRef = useRef<Map<string, SimNode>>(new Map());
  const isSettledRef = useRef(false);
  const iterationRef = useRef(0);
  const alphaRef = useRef(FORCES.initialAlpha);

  // Convert elements to simulation nodes/edges
  useMemo(() => {
    const newNodes: SimNode[] = [];
    const newEdges: SimEdge[] = [];
    const nodeMap = new Map<string, SimNode>();

    // Get existing positions if available
    const oldNodeMap = nodeMapRef.current;

    // Package order for consistent angular placement
    const packageOrder: PackageId[] = [
      "numwasm",
      "sciwasm",
      "symwasm",
      "arwasm",
      "lawasm",
      "linwasm",
      "quadwasm",
      "superluwasm",
      "xsfwasm",
    ];

    // Generate package centers dynamically at fixed distance from origin (spherical distribution)
    const packageCenters: Partial<
      Record<PackageId, { x: number; y: number; z: number }>
    > = {};
    for (let i = 0; i < packageOrder.length; i++) {
      // Distribute on a sphere using golden angle for even spacing
      const phi = Math.acos(1 - (2 * (i + 0.5)) / packageOrder.length);
      const theta = Math.PI * (1 + Math.sqrt(5)) * i;
      packageCenters[packageOrder[i]] = {
        x: Math.sin(phi) * Math.cos(theta) * FORCES.packageDistance,
        y: Math.sin(phi) * Math.sin(theta) * FORCES.packageDistance,
        z: Math.cos(phi) * FORCES.packageDistance,
      };
    }

    // Separate elements by type
    const packageElements: (CyElement & { data: CyNodeData })[] = [];
    const wasmModuleElements: (CyElement & { data: CyNodeData })[] = [];
    const functionElements: (CyElement & { data: CyNodeData })[] = [];
    const edgeElements: (CyElement & { data: CyEdgeData })[] = [];

    for (const el of elements) {
      if (isNodeElement(el)) {
        if (el.data.type === "package") packageElements.push(el);
        else if (el.data.type === "wasmModule") wasmModuleElements.push(el);
        else functionElements.push(el);
      } else if (isEdgeElement(el)) {
        edgeElements.push(el);
      }
    }

    // Helper to create a node
    const createNode = (
      data: CyNodeData,
      x: number,
      y: number,
      z: number,
    ): SimNode => ({
      id: data.id,
      x,
      y,
      z,
      vx: 0,
      vy: 0,
      vz: 0,
      type: data.type,
      packageId: data.packageId as PackageId,
      color: data.color,
      label: data.label,
      mass: NODE_MASS[data.type] || 1,
      radius: NODE_RADIUS[data.type] || 0.5,
      kind: data.kind,
      isShared: data.isShared,
      description: data.description,
      searchScore: 0, // Updated separately by searchQuery effect
    });

    // Phase 1: Place packages at their centers
    for (const el of packageElements) {
      const data = el.data;
      const oldNode = oldNodeMap.get(data.id);
      const center = packageCenters[data.packageId as PackageId] || {
        x: 0,
        y: 0,
        z: 0,
      };

      const node = createNode(
        data,
        oldNode?.x ?? center.x,
        oldNode?.y ?? center.y,
        oldNode?.z ?? center.z,
      );
      newNodes.push(node);
      nodeMap.set(data.id, node);
    }

    // Phase 2: Place WASM modules in a ring around their package
    for (let i = 0; i < wasmModuleElements.length; i++) {
      const el = wasmModuleElements[i];
      const data = el.data;
      const oldNode = oldNodeMap.get(data.id);
      const packageNode = nodeMap.get(`pkg-${data.packageId}`);
      const center = packageNode ||
        packageCenters[data.packageId as PackageId] || { x: 0, y: 0, z: 0 };

      // Count modules per package to distribute evenly
      const samePackageModules = wasmModuleElements.filter(
        (e) => e.data.packageId === data.packageId,
      );
      const indexInPackage = samePackageModules.indexOf(el);
      const countInPackage = samePackageModules.length;

      // Distribute in a circle around the package
      const angle = (indexInPackage / countInPackage) * Math.PI * 2;
      const radius = 20; // Distance from package center

      const node = createNode(
        data,
        oldNode?.x ?? center.x + Math.cos(angle) * radius,
        oldNode?.y ?? center.y + Math.sin(angle) * radius,
        oldNode?.z ?? center.z + (Math.random() - 0.5) * 2,
      );
      newNodes.push(node);
      nodeMap.set(data.id, node);
    }

    // Add all edges
    for (const el of edgeElements) {
      const data = el.data;
      newEdges.push({
        source: data.source,
        target: data.target,
        edgeType: data.edgeType,
        color: data.color,
      });
    }

    // Run a few iterations with just packages and WASM modules
    const WARMUP_ITERATIONS = 15;
    const tempNodes = newNodes;
    const tempEdges = newEdges.filter((e) => {
      const source = nodeMap.get(e.source);
      const target = nodeMap.get(e.target);
      return (
        source &&
        target &&
        source.type !== "function" &&
        target.type !== "function"
      );
    });

    for (let iter = 0; iter < WARMUP_ITERATIONS; iter++) {
      const alpha = 1.0 - (iter / WARMUP_ITERATIONS) * 0.5;

      // Apply forces to packages and WASM modules only
      for (const node of tempNodes) {
        let fx = 0,
          fy = 0,
          fz = 0;

        // Repulsion
        for (const other of tempNodes) {
          if (node.id === other.id) continue;
          const dx = node.x - other.x;
          const dy = node.y - other.y;
          const dz = node.z - other.z;
          const distSq = dx * dx + dy * dy + dz * dz;
          const dist = Math.sqrt(distSq) || 0.1;
          const force = (FORCES.repulsion * 2 * alpha) / (distSq + 1);
          fx += (dx / dist) * force;
          fy += (dy / dist) * force;
          fz += (dz / dist) * force;
        }

        // Center pull
        fx -= node.x * FORCES.centerPull * alpha;
        fy -= node.y * FORCES.centerPull * alpha;
        fz -= node.z * FORCES.centerPull * alpha;

        node.vx += fx / node.mass;
        node.vy += fy / node.mass;
        node.vz += fz / node.mass;
      }

      // Edge attraction
      for (const edge of tempEdges) {
        const source = nodeMap.get(edge.source);
        const target = nodeMap.get(edge.target);
        if (!source || !target) continue;

        const dx = target.x - source.x;
        const dy = target.y - source.y;
        const dz = target.z - source.z;
        const dist = Math.sqrt(dx * dx + dy * dy + dz * dz) || 0.1;

        const idealDist = source.type === "package" ? 25 : 15;
        const displacement = dist - idealDist;
        const force = displacement * FORCES.attraction * alpha;

        source.vx += ((dx / dist) * force) / source.mass;
        source.vy += ((dy / dist) * force) / source.mass;
        target.vx -= ((dx / dist) * force) / target.mass;
        target.vy -= ((dy / dist) * force) / target.mass;
      }

      // Update positions
      for (const node of tempNodes) {
        node.x += node.vx;
        node.y += node.vy;
        node.z += node.vz;
        node.vx *= 0.8;
        node.vy *= 0.8;
        node.vz *= 0.8;
      }
    }

    // Phase 3: Now place function nodes near their parent (WASM module or package)
    for (const el of functionElements) {
      const data = el.data;
      const oldNode = oldNodeMap.get(data.id);

      // Find the parent node (WASM module for wasm functions, package for TS functions)
      let parentNode: SimNode | undefined;
      if (data.kind === "wasm") {
        // Find the WASM module this function belongs to
        const moduleEdge = newEdges.find(
          (e) => e.target === data.id && e.edgeType === "containment",
        );
        if (moduleEdge) {
          parentNode = nodeMap.get(moduleEdge.source);
        }
      }
      // Fall back to package node
      if (!parentNode) {
        parentNode = nodeMap.get(`pkg-${data.packageId}`);
      }

      const center = parentNode ||
        packageCenters[data.packageId as PackageId] || { x: 0, y: 0, z: 0 };

      // Spawn in a small neighborhood around the parent
      const theta = Math.random() * Math.PI * 2;
      const phi = Math.acos(2 * Math.random() - 1);
      const r = 8 * Math.cbrt(Math.random());

      const node = createNode(
        data,
        oldNode?.x ?? center.x + r * Math.sin(phi) * Math.cos(theta),
        oldNode?.y ?? center.y + r * Math.sin(phi) * Math.sin(theta),
        oldNode?.z ?? (center.z || 0) + r * Math.cos(phi) * 0.2,
      );
      newNodes.push(node);
      nodeMap.set(data.id, node);
    }

    nodesRef.current = newNodes;
    edgesRef.current = newEdges;
    nodeMapRef.current = nodeMap;
    isSettledRef.current = false;
    iterationRef.current = 0;
    alphaRef.current = FORCES.initialAlpha;
  }, [elements]); // Note: searchQuery intentionally excluded - handled separately

  // Update search scores without resetting positions
  useMemo(() => {
    const nodes = nodesRef.current;
    for (const node of nodes) {
      node.searchScore = searchQuery
        ? calculateSearchScore(searchQuery, node.label)
        : 0;
    }
  }, [searchQuery]);

  // Force simulation step
  const step = useCallback(() => {
    const nodes = nodesRef.current;
    const edges = edgesRef.current;
    const nodeMap = nodeMapRef.current;

    if (nodes.length === 0) return;

    // Get current alpha (temperature)
    const alpha = alphaRef.current;

    // Stop when alpha is too low (simulation has cooled)
    if (alpha < FORCES.alphaMin) {
      isSettledRef.current = true;
      return;
    }

    // Decay alpha (simulated annealing)
    alphaRef.current = alpha * (1 - FORCES.alphaDecay);
    iterationRef.current++;

    // Pre-compute node counts per package for normalized cluster pull
    const packageNodeCounts = new Map<PackageId, number>();
    for (const node of nodes) {
      if (node.type !== "package") {
        const count = packageNodeCounts.get(node.packageId) || 0;
        packageNodeCounts.set(node.packageId, count + 1);
      }
    }

    // Calculate forces (scaled by alpha for annealing)
    for (const node of nodes) {
      let fx = 0,
        fy = 0,
        fz = 0;

      // 1. Repulsion from all other nodes (optimized: only check nearby)
      for (const other of nodes) {
        if (node.id === other.id) continue;

        const dx = node.x - other.x;
        const dy = node.y - other.y;
        const dz = node.z - other.z;
        const distSq = dx * dx + dy * dy + dz * dz;
        const dist = Math.sqrt(distSq) || FORCES.minDistance;

        // Strong repulsion between packages (no distance limit)
        if (node.type === "package" && other.type === "package") {
          const force = (FORCES.repulsion * 50 * alpha) / (distSq + 1);
          fx += (dx / dist) * force;
          fy += (dy / dist) * force;
          fz += (dz / dist) * force;
        } else if (dist < 50) {
          // Only compute for nearby nodes
          let force = (FORCES.repulsion * alpha) / (distSq + 1);
          fx += (dx / dist) * force;
          fy += (dy / dist) * force;
          fz += (dz / dist) * force;

          if (node.kind === "wasm" && other.kind === "wasm") {
            force *= 20;
          }
        }
      }

      // 2. Center pull (gentle gravity toward origin)
      fx -= node.x * FORCES.centerPull * alpha;
      fy -= node.y * FORCES.centerPull * alpha;
      fz -= node.z * FORCES.centerPull * alpha;

      // 3. Cluster pull (nodes with same packageId attract, normalized by package size)
      const packageCount = packageNodeCounts.get(node.packageId) || 1;
      for (const other of nodes) {
        if (node.id === other.id) continue;
        if (node.packageId === other.packageId && node.type !== "package") {
          const dx = other.x - node.x;
          const dy = other.y - node.y;
          const dz = other.z - node.z;
          const dist =
            Math.sqrt(dx * dx + dy * dy + dz * dz) || FORCES.minDistance;

          // Normalize by package node count so total attraction is independent of package size
          const normalizedPull =
            (FORCES.clusterPull * alpha * dist * 0.05) / packageCount;
          fx += (dx / dist) * normalizedPull;
          fy += (dy / dist) * normalizedPull;
          fz += (dz / dist) * normalizedPull;
        }
      }

      // Apply forces (F = ma, so a = F/m)
      node.vx += fx / node.mass;
      node.vy += fy / node.mass;
      node.vz += fz / node.mass;
    }

    // 4. Edge attraction (spring force, also scaled by alpha)
    for (const edge of edges) {
      const source = nodeMap.get(edge.source);
      const target = nodeMap.get(edge.target);
      if (!source || !target) continue;

      const dx = target.x - source.x;
      const dy = target.y - source.y;
      const dz = target.z - source.z;
      const dist = Math.sqrt(dx * dx + dy * dy + dz * dz) || FORCES.minDistance;

      // Determine attraction strength based on edge and node types
      let attractionStrength = FORCES.attraction;
      let idealDist = (source.radius + target.radius) * 4;

      if (edge.edgeType === "binding") {
        attractionStrength = 2 * FORCES.bindingAttraction;
        idealDist = (source.radius + target.radius) * 2; // Closer together
      } else if (source.type === "package" && target.type === "wasmModule") {
        // Package to WASM module: prefer long distance, acts more like repulsion
        attractionStrength = FORCES.attraction * 0.001;
        idealDist = 100; // Push modules far from package center
      } else if (source.type === "package" || target.type === "package") {
        // Other package edges (package to TS functions)
        attractionStrength = FORCES.attraction * 0.5;
        idealDist = 20;
      }

      const displacement = dist - idealDist;
      const force = displacement * attractionStrength * alpha;

      const fx = (dx / dist) * force;
      const fy = (dy / dist) * force;
      const fz = (dz / dist) * force;

      source.vx += fx / source.mass;
      source.vy += fy / source.mass;
      source.vz += fz / source.mass;
      target.vx -= fx / target.mass;
      target.vy -= fy / target.mass;
      target.vz -= fz / target.mass;
    }

    // Apply velocities and damping
    for (const node of nodes) {
      // Clamp velocity (scaled by alpha for smoother settling)
      const maxVel = FORCES.maxVelocity * Math.max(alpha, 0.1);
      const speed = Math.sqrt(
        node.vx * node.vx + node.vy * node.vy + node.vz * node.vz,
      );
      if (speed > maxVel) {
        const scale = maxVel / speed;
        node.vx *= scale;
        node.vy *= scale;
        node.vz *= scale;
      }

      // Update position
      node.x += node.vx;
      node.y += node.vy;
      node.z += node.vz;

      // Hard constraint: packages must stay at fixed distance from center
      if (node.type === "package") {
        const currentDist =
          Math.sqrt(node.x * node.x + node.y * node.y + node.z * node.z) ||
          FORCES.minDistance;
        const scale = FORCES.packageDistance / currentDist;
        node.x *= scale;
        node.y *= scale;
        node.z *= scale;
      }

      // Apply velocity decay (in addition to damping)
      const decay = FORCES.damping * (1 - FORCES.velocityDecay * (1 - alpha));
      node.vx *= decay;
      node.vy *= decay;
      node.vz *= decay;
    }
  }, []);

  return {
    nodes: nodesRef.current,
    edges: edgesRef.current,
    step,
    isSettled: isSettledRef.current,
  };
}
