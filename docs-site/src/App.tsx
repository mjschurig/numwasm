import { useState } from "react";
import { Routes, Route, NavLink } from "react-router-dom";
import { Link, Button } from "react-aria-components";
import { Menu, X } from "lucide-react";
import Home from "./pages/Home";
import Demo from "./pages/Demo";
import Docs from "./pages/Docs";
import Benchmarks from "./pages/Benchmarks";
import { Background } from "./components/background";

// GitHub icon from Simple Icons (https://simpleicons.org)
function GitHubIcon({ className }: { className?: string }) {
  return (
    <svg
      role="img"
      viewBox="0 0 24 24"
      fill="currentColor"
      className={className}
      aria-hidden="true"
    >
      <path d="M12 .297c-6.63 0-12 5.373-12 12 0 5.303 3.438 9.8 8.205 11.385.6.113.82-.258.82-.577 0-.285-.01-1.04-.015-2.04-3.338.724-4.042-1.61-4.042-1.61C4.422 18.07 3.633 17.7 3.633 17.7c-1.087-.744.084-.729.084-.729 1.205.084 1.838 1.236 1.838 1.236 1.07 1.835 2.809 1.305 3.495.998.108-.776.417-1.305.76-1.605-2.665-.3-5.466-1.332-5.466-5.93 0-1.31.465-2.38 1.235-3.22-.135-.303-.54-1.523.105-3.176 0 0 1.005-.322 3.3 1.23.96-.267 1.98-.399 3-.405 1.02.006 2.04.138 3 .405 2.28-1.552 3.285-1.23 3.285-1.23.645 1.653.24 2.873.12 3.176.765.84 1.23 1.91 1.23 3.22 0 4.61-2.805 5.625-5.475 5.92.42.36.81 1.096.81 2.22 0 1.606-.015 2.896-.015 3.286 0 .315.21.69.825.57C20.565 22.092 24 17.592 24 12.297c0-6.627-5.373-12-12-12" />
    </svg>
  );
}

function Navigation() {
  const [mobileMenuOpen, setMobileMenuOpen] = useState(false);

  const navLinkClass = ({ isActive }: { isActive: boolean }) =>
    `font-medium py-2 border-b-2 transition-colors no-underline hover:no-underline ${
      isActive
        ? "border-primary text-primary"
        : "border-transparent text-gray-300 hover:text-primary"
    }`;

  const mobileNavLinkClass = ({ isActive }: { isActive: boolean }) =>
    `block px-4 py-3 font-medium transition-colors no-underline hover:no-underline ${
      isActive
        ? "text-primary bg-gray-800/50"
        : "text-gray-300 hover:text-primary hover:bg-gray-800/30"
    }`;

  return (
    <nav className="bg-[#030712]/80 backdrop-blur-sm border-b border-gray-800 sticky top-0 z-50">
      <div className="flex items-center justify-between px-4 sm:px-8 h-16">
        <NavLink
          to="/"
          className="text-xl sm:text-2xl font-bold hero-text no-underline hover:no-underline"
        >
          numwasm
        </NavLink>

        {/* Desktop nav */}
        <div className="hidden md:flex items-center gap-8">
          <NavLink to="/" end className={navLinkClass}>
            Home
          </NavLink>
          <NavLink to="/demo" className={navLinkClass}>
            Demo
          </NavLink>
          <NavLink to="/docs" className={navLinkClass}>
            Docs
          </NavLink>
          <NavLink to="/benchmarks" className={navLinkClass}>
            Benchmarks
          </NavLink>
          <Link
            href="https://github.com/mjschurig/numwasm"
            target="_blank"
            rel="noopener noreferrer"
            className="flex items-center gap-2 font-medium text-gray-300 hover:text-primary no-underline"
          >
            <GitHubIcon className="w-5 h-5" />
            GitHub
          </Link>
        </div>

        {/* Mobile hamburger button */}
        <Button
          onPress={() => setMobileMenuOpen(!mobileMenuOpen)}
          className="md:hidden p-2 text-gray-300 hover:text-primary cursor-pointer"
          aria-label="Toggle menu"
        >
          {mobileMenuOpen ? (
            <X className="w-6 h-6" />
          ) : (
            <Menu className="w-6 h-6" />
          )}
        </Button>
      </div>

      {/* Mobile menu */}
      {mobileMenuOpen && (
        <div className="md:hidden border-t border-gray-800 bg-[#030712]/95 backdrop-blur-sm">
          <NavLink
            to="/"
            end
            className={mobileNavLinkClass}
            onClick={() => setMobileMenuOpen(false)}
          >
            Home
          </NavLink>
          <NavLink
            to="/demo"
            className={mobileNavLinkClass}
            onClick={() => setMobileMenuOpen(false)}
          >
            Demo
          </NavLink>
          <NavLink
            to="/docs"
            className={mobileNavLinkClass}
            onClick={() => setMobileMenuOpen(false)}
          >
            Docs
          </NavLink>
          <NavLink
            to="/benchmarks"
            className={mobileNavLinkClass}
            onClick={() => setMobileMenuOpen(false)}
          >
            Benchmarks
          </NavLink>
          <Link
            href="https://github.com/mjschurig/numwasm"
            target="_blank"
            rel="noopener noreferrer"
            className="flex items-center gap-2 px-4 py-3 font-medium text-gray-300 hover:text-primary hover:bg-gray-800/30 no-underline"
          >
            <GitHubIcon className="w-5 h-5" />
            GitHub
          </Link>
        </div>
      )}
    </nav>
  );
}

export default function App() {
  return (
    <div className="relative min-h-screen flex flex-col">
      <div className="fixed inset-0 z-0 overflow-hidden">
        <Background />
      </div>
      <div className="relative z-10 flex flex-col min-h-screen">
        <Navigation />
        <div className="flex-1">
          <Routes>
            <Route path="/" element={<Home />} />
            <Route path="/demo" element={<Demo />} />
            <Route path="/docs/*" element={<Docs />} />
            <Route path="/benchmarks/*" element={<Benchmarks />} />
          </Routes>
        </div>
      </div>
    </div>
  );
}
