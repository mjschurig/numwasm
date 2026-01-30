import { useRef, useCallback, useMemo } from 'react';
import type { CyElement, CyNodeData, CyEdgeData, PackageId } from './types';

export interface SimNode {
  id: string;
  x: number;
  y: number;
  z: number;
  vx: number;
  vy: number;
  vz: number;
  type: 'package' | 'wasmModule' | 'function';
  packageId: PackageId;
  color: string;
  label: string;
  mass: number;
  radius: number;
  kind?: 'wasm' | 'typescript';
  isShared?: boolean;
}

export interface SimEdge {
  source: string;
  target: string;
  edgeType: 'containment' | 'binding';
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
  repulsion: 800,
  attraction: 0.03,
  centerPull: 0.005,
  clusterPull: 0.01,
  damping: 0.85,
  maxVelocity: 5,
  minDistance: 0.5,
};

// Node sizes by type
const NODE_RADIUS: Record<string, number> = {
  package: 2.0,
  wasmModule: 1.2,
  function: 0.5,
};

const NODE_MASS: Record<string, number> = {
  package: 10,
  wasmModule: 5,
  function: 1,
};

function isNodeElement(el: CyElement): el is CyElement & { data: CyNodeData } {
  return !('source' in el.data);
}

function isEdgeElement(el: CyElement): el is CyElement & { data: CyEdgeData } {
  return 'source' in el.data;
}

export function useForceSimulation(elements: CyElement[]): ForceSimulationResult {
  const nodesRef = useRef<SimNode[]>([]);
  const edgesRef = useRef<SimEdge[]>([]);
  const nodeMapRef = useRef<Map<string, SimNode>>(new Map());
  const isSettledRef = useRef(false);
  const iterationRef = useRef(0);

  // Convert elements to simulation nodes/edges
  useMemo(() => {
    const newNodes: SimNode[] = [];
    const newEdges: SimEdge[] = [];
    const nodeMap = new Map<string, SimNode>();

    // Get existing positions if available
    const oldNodeMap = nodeMapRef.current;

    for (const el of elements) {
      if (isNodeElement(el)) {
        const data = el.data;
        const oldNode = oldNodeMap.get(data.id);

        // Initial position: use old position or random
        const spread = data.type === 'package' ? 20 : data.type === 'wasmModule' ? 40 : 60;
        const x = oldNode?.x ?? (Math.random() - 0.5) * spread;
        const y = oldNode?.y ?? (Math.random() - 0.5) * spread;
        const z = oldNode?.z ?? (Math.random() - 0.5) * spread * 0.3;

        const node: SimNode = {
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
        };
        newNodes.push(node);
        nodeMap.set(data.id, node);
      } else if (isEdgeElement(el)) {
        const data = el.data;
        newEdges.push({
          source: data.source,
          target: data.target,
          edgeType: data.edgeType,
          color: data.color,
        });
      }
    }

    nodesRef.current = newNodes;
    edgesRef.current = newEdges;
    nodeMapRef.current = nodeMap;
    isSettledRef.current = false;
    iterationRef.current = 0;
  }, [elements]);

  // Force simulation step
  const step = useCallback(() => {
    const nodes = nodesRef.current;
    const edges = edgesRef.current;
    const nodeMap = nodeMapRef.current;

    if (nodes.length === 0) return;

    // Stop after enough iterations
    if (iterationRef.current > 500) {
      isSettledRef.current = true;
      return;
    }
    iterationRef.current++;

    // Calculate forces
    for (const node of nodes) {
      let fx = 0, fy = 0, fz = 0;

      // 1. Repulsion from all other nodes (optimized: only check nearby)
      for (const other of nodes) {
        if (node.id === other.id) continue;

        const dx = node.x - other.x;
        const dy = node.y - other.y;
        const dz = node.z - other.z;
        const distSq = dx * dx + dy * dy + dz * dz;
        const dist = Math.sqrt(distSq) || FORCES.minDistance;

        if (dist < 50) { // Only compute for nearby nodes
          const force = FORCES.repulsion / (distSq + 1);
          fx += (dx / dist) * force;
          fy += (dy / dist) * force;
          fz += (dz / dist) * force * 0.3; // Less force in Z
        }
      }

      // 2. Center pull (gentle gravity toward origin)
      fx -= node.x * FORCES.centerPull;
      fy -= node.y * FORCES.centerPull;
      fz -= node.z * FORCES.centerPull * 2; // Stronger in Z to keep flat

      // 3. Cluster pull (nodes with same packageId attract)
      for (const other of nodes) {
        if (node.id === other.id) continue;
        if (node.packageId === other.packageId && node.type !== 'package') {
          const dx = other.x - node.x;
          const dy = other.y - node.y;
          const dz = other.z - node.z;
          const dist = Math.sqrt(dx * dx + dy * dy + dz * dz) || FORCES.minDistance;

          fx += (dx / dist) * FORCES.clusterPull * dist * 0.1;
          fy += (dy / dist) * FORCES.clusterPull * dist * 0.1;
          fz += (dz / dist) * FORCES.clusterPull * dist * 0.05;
        }
      }

      // Apply forces (F = ma, so a = F/m)
      node.vx += fx / node.mass;
      node.vy += fy / node.mass;
      node.vz += fz / node.mass;
    }

    // 4. Edge attraction (spring force)
    for (const edge of edges) {
      const source = nodeMap.get(edge.source);
      const target = nodeMap.get(edge.target);
      if (!source || !target) continue;

      const dx = target.x - source.x;
      const dy = target.y - source.y;
      const dz = target.z - source.z;
      const dist = Math.sqrt(dx * dx + dy * dy + dz * dz) || FORCES.minDistance;

      // Ideal distance based on node types
      const idealDist = (source.radius + target.radius) * 4;
      const displacement = dist - idealDist;
      const force = displacement * FORCES.attraction;

      const fx = (dx / dist) * force;
      const fy = (dy / dist) * force;
      const fz = (dz / dist) * force * 0.3;

      source.vx += fx / source.mass;
      source.vy += fy / source.mass;
      source.vz += fz / source.mass;
      target.vx -= fx / target.mass;
      target.vy -= fy / target.mass;
      target.vz -= fz / target.mass;
    }

    // Apply velocities and damping
    let totalMovement = 0;
    for (const node of nodes) {
      // Clamp velocity
      const speed = Math.sqrt(node.vx * node.vx + node.vy * node.vy + node.vz * node.vz);
      if (speed > FORCES.maxVelocity) {
        const scale = FORCES.maxVelocity / speed;
        node.vx *= scale;
        node.vy *= scale;
        node.vz *= scale;
      }

      // Update position
      node.x += node.vx;
      node.y += node.vy;
      node.z += node.vz;

      // Damping
      node.vx *= FORCES.damping;
      node.vy *= FORCES.damping;
      node.vz *= FORCES.damping;

      totalMovement += Math.abs(node.vx) + Math.abs(node.vy) + Math.abs(node.vz);
    }

    // Check if settled
    if (totalMovement < 0.1 * nodes.length) {
      isSettledRef.current = true;
    }
  }, []);

  return {
    nodes: nodesRef.current,
    edges: edgesRef.current,
    step,
    isSettled: isSettledRef.current,
  };
}
