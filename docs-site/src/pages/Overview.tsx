import SEO from '../components/SEO';
import { OverviewGraph } from '../components/overview';

export default function Overview() {
  return (
    <>
      <SEO
        title="Overview"
        description="Interactive visualization of WASM modules and TypeScript exports across numwasm, sciwasm, and symwasm packages."
        path="/overview"
      />
      <div className="h-[calc(100vh-4rem)] flex flex-col">
        <OverviewGraph />
      </div>
    </>
  );
}
