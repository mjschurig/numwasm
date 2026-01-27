import MatrixDemo from '../components/demo/MatrixDemo';
import SEO from '../components/SEO';

export default function Demo() {
  return (
    <>
      <SEO
        title="Interactive Demo"
        description="Try numwasm in the browser. Interactive matrix operations demo with WebAssembly acceleration."
        path="/demo"
        breadcrumbs={[
          { name: 'Home', url: '/' },
          { name: 'Demo', url: '/demo' },
        ]}
      />
      <MatrixDemo />
    </>
  );
}
