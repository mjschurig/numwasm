// TypeDoc JSON output types
// Based on TypeDoc's JSONOutput module

export interface ProjectReflection {
  id: number;
  name: string;
  kind: number;
  children?: DeclarationReflection[];
  groups?: ReflectionGroup[];
  symbolIdMap?: Record<string, number>;
}

export interface ReflectionGroup {
  title: string;
  children: number[];
}

export interface DeclarationReflection {
  id: number;
  name: string;
  kind: number;
  kindString?: string;
  flags?: Flags;
  comment?: Comment;
  children?: DeclarationReflection[];
  groups?: ReflectionGroup[];
  signatures?: SignatureReflection[];
  type?: TypeInfo;
  defaultValue?: string;
  sources?: SourceReference[];
  extendedTypes?: TypeInfo[];
  implementedTypes?: TypeInfo[];
  typeParameters?: TypeParameterReflection[];
  getSignature?: SignatureReflection;
  setSignature?: SignatureReflection;
}

export interface SignatureReflection {
  id: number;
  name: string;
  kind: number;
  kindString?: string;
  flags?: Flags;
  comment?: Comment;
  parameters?: ParameterReflection[];
  type?: TypeInfo;
  typeParameter?: TypeParameterReflection[];
}

export interface ParameterReflection {
  id: number;
  name: string;
  kind: number;
  kindString?: string;
  flags?: Flags;
  comment?: Comment;
  type?: TypeInfo;
  defaultValue?: string;
}

export interface TypeParameterReflection {
  id: number;
  name: string;
  kind: number;
  type?: TypeInfo;
  default?: TypeInfo;
}

export interface Flags {
  isPrivate?: boolean;
  isProtected?: boolean;
  isPublic?: boolean;
  isStatic?: boolean;
  isExported?: boolean;
  isExternal?: boolean;
  isOptional?: boolean;
  isRest?: boolean;
  isAbstract?: boolean;
  isConst?: boolean;
  isReadonly?: boolean;
}

export interface Comment {
  summary?: CommentDisplayPart[];
  blockTags?: CommentTag[];
  modifierTags?: string[];
}

export interface CommentDisplayPart {
  kind: 'text' | 'code' | 'inline-tag';
  text: string;
  tag?: string;
  target?: number | string;
}

export interface CommentTag {
  tag: string;
  content: CommentDisplayPart[];
}

export interface TypeInfo {
  type: string;
  name?: string;
  value?: string | number | boolean | null;
  elementType?: TypeInfo;
  types?: TypeInfo[];
  typeArguments?: TypeInfo[];
  declaration?: DeclarationReflection;
  target?: number;
  id?: number;
  qualifiedName?: string;
  package?: string;
}

export interface SourceReference {
  fileName: string;
  line: number;
  character: number;
  url?: string;
}

// Kind enum values from TypeDoc
export const ReflectionKind = {
  Project: 1,
  Module: 2,
  Namespace: 4,
  Enum: 8,
  EnumMember: 16,
  Variable: 32,
  Function: 64,
  Class: 128,
  Interface: 256,
  Constructor: 512,
  Property: 1024,
  Method: 2048,
  CallSignature: 4096,
  IndexSignature: 8192,
  ConstructorSignature: 16384,
  Parameter: 32768,
  TypeLiteral: 65536,
  TypeParameter: 131072,
  Accessor: 262144,
  GetSignature: 524288,
  SetSignature: 1048576,
  TypeAlias: 2097152,
  Reference: 4194304,
} as const;

export function getKindString(kind: number): string {
  switch (kind) {
    case ReflectionKind.Module: return 'Module';
    case ReflectionKind.Namespace: return 'Namespace';
    case ReflectionKind.Enum: return 'Enum';
    case ReflectionKind.EnumMember: return 'Enum Member';
    case ReflectionKind.Variable: return 'Variable';
    case ReflectionKind.Function: return 'Function';
    case ReflectionKind.Class: return 'Class';
    case ReflectionKind.Interface: return 'Interface';
    case ReflectionKind.Constructor: return 'Constructor';
    case ReflectionKind.Property: return 'Property';
    case ReflectionKind.Method: return 'Method';
    case ReflectionKind.TypeAlias: return 'Type Alias';
    case ReflectionKind.Accessor: return 'Accessor';
    default: return 'Unknown';
  }
}
