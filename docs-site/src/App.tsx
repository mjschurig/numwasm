import { Routes, Route, NavLink } from 'react-router-dom';
import { Link } from 'react-aria-components';
import Home from './pages/Home';
import Demo from './pages/Demo';
import Docs from './pages/Docs';
import { Background } from './components/background';

function Navigation() {
  return (
    <nav className="flex items-center justify-between px-8 h-16 bg-[#030712]/80 backdrop-blur-sm border-b border-gray-800 sticky top-0 z-50">
      <NavLink to="/" className="text-2xl font-bold hero-text no-underline hover:no-underline">
        numwasm
      </NavLink>
      <div className="flex items-center gap-8">
        <NavLink
          to="/"
          end
          className={({ isActive }) =>
            `font-medium py-2 border-b-2 transition-colors no-underline hover:no-underline ${
              isActive ? 'border-primary text-primary' : 'border-transparent text-gray-300 hover:text-primary'
            }`
          }
        >
          Home
        </NavLink>
        <NavLink
          to="/demo"
          className={({ isActive }) =>
            `font-medium py-2 border-b-2 transition-colors no-underline hover:no-underline ${
              isActive ? 'border-primary text-primary' : 'border-transparent text-gray-300 hover:text-primary'
            }`
          }
        >
          Demo
        </NavLink>
        <NavLink
          to="/docs"
          className={({ isActive }) =>
            `font-medium py-2 border-b-2 transition-colors no-underline hover:no-underline ${
              isActive ? 'border-primary text-primary' : 'border-transparent text-gray-300 hover:text-primary'
            }`
          }
        >
          Docs
        </NavLink>
        <Link
          href="https://github.com/nicholasgriffintn/numwasm"
          target="_blank"
          rel="noopener noreferrer"
          className="font-medium text-gray-300 hover:text-primary no-underline"
        >
          GitHub
        </Link>
      </div>
    </nav>
  );
}

export default function App() {
  return (
    <div className="relative min-h-screen overflow-hidden">
      <Background />
      <div className="relative z-10">
        <Navigation />
        <Routes>
          <Route path="/" element={<Home />} />
          <Route path="/demo" element={<Demo />} />
          <Route path="/docs/*" element={<Docs />} />
        </Routes>
      </div>
    </div>
  );
}
