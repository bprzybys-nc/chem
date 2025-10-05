import Link from 'next/link'

export default function Home() {
  return (
    <div className="min-h-screen flex items-center justify-center bg-gradient-to-br from-blue-50 to-indigo-100">
      <div className="text-center">
        <h1 className="text-6xl font-bold mb-4 text-gray-900">
          Chemistry 3D Visualizer
        </h1>
        <p className="text-xl text-gray-600 mb-8">
          Interactive molecular visualization powered by 3Dmol.js and RDKit
        </p>
        <div className="flex gap-4 justify-center">
          <Link
            href="/visualize"
            className="px-8 py-3 bg-blue-600 text-white rounded-lg hover:bg-blue-700 transition-colors text-lg font-semibold shadow-lg"
          >
            Start Visualizing
          </Link>
          <a
            href="http://localhost:8000/docs"
            target="_blank"
            rel="noopener noreferrer"
            className="px-8 py-3 bg-white text-blue-600 border-2 border-blue-600 rounded-lg hover:bg-blue-50 transition-colors text-lg font-semibold shadow-lg"
          >
            API Documentation
          </a>
        </div>
        <div className="mt-12 text-gray-500">
          <p className="text-sm">Built with Next.js, FastAPI, and RDKit</p>
        </div>
      </div>
    </div>
  )
}
