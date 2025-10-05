'use client'

import { useState } from 'react'
import MoleculeViewer from '@/components/MoleculeViewer'

export default function VisualizePage() {
  const [smiles, setSmiles] = useState('c1ccccc1') // Benzene as default
  const [pdbData, setPdbData] = useState('')
  const [loading, setLoading] = useState(false)
  const [error, setError] = useState('')

  const handleVisualize = async () => {
    if (!smiles.trim()) {
      setError('Please enter a SMILES string')
      return
    }

    setLoading(true)
    setError('')

    try {
      // In Phase 3, this will call the backend API
      // For now, we'll use a hardcoded benzene PDB
      const benzenePDB = `COMPND    benzene
HETATM    1  C   UNL     1       1.200   0.000   0.000  1.00  0.00           C
HETATM    2  C   UNL     1       0.600   1.039   0.000  1.00  0.00           C
HETATM    3  C   UNL     1      -0.600   1.039   0.000  1.00  0.00           C
HETATM    4  C   UNL     1      -1.200   0.000   0.000  1.00  0.00           C
HETATM    5  C   UNL     1      -0.600  -1.039   0.000  1.00  0.00           C
HETATM    6  C   UNL     1       0.600  -1.039   0.000  1.00  0.00           C
HETATM    7  H   UNL     1       2.130   0.000   0.000  1.00  0.00           H
HETATM    8  H   UNL     1       1.065   1.845   0.000  1.00  0.00           H
HETATM    9  H   UNL     1      -1.065   1.845   0.000  1.00  0.00           H
HETATM   10  H   UNL     1      -2.130   0.000   0.000  1.00  0.00           H
HETATM   11  H   UNL     1      -1.065  -1.845   0.000  1.00  0.00           H
HETATM   12  H   UNL     1       1.065  -1.845   0.000  1.00  0.00           H
CONECT    1    2    6    7
CONECT    2    1    3    8
CONECT    3    2    4    9
CONECT    4    3    5   10
CONECT    5    4    6   11
CONECT    6    1    5   12
END`

      setPdbData(benzenePDB)
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to visualize molecule')
    } finally {
      setLoading(false)
    }
  }

  return (
    <div className="min-h-screen p-8 bg-gray-50">
      <div className="max-w-6xl mx-auto">
        <h1 className="text-4xl font-bold mb-8 text-gray-900">
          3D Molecule Visualizer
        </h1>

        <div className="bg-white p-6 rounded-lg shadow-md mb-8">
          <label htmlFor="smiles" className="block text-sm font-medium text-gray-700 mb-2">
            SMILES String
          </label>
          <div className="flex gap-4">
            <input
              id="smiles"
              type="text"
              value={smiles}
              onChange={(e) => setSmiles(e.target.value)}
              placeholder="Enter SMILES (e.g., c1ccccc1 for benzene)"
              className="flex-1 px-4 py-2 border border-gray-300 rounded-md focus:ring-2 focus:ring-blue-500 focus:border-blue-500"
            />
            <button
              onClick={handleVisualize}
              disabled={loading}
              className="px-6 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700 disabled:bg-gray-400 disabled:cursor-not-allowed transition-colors"
            >
              {loading ? 'Loading...' : 'Visualize'}
            </button>
          </div>
          {error && (
            <p className="mt-2 text-sm text-red-600">{error}</p>
          )}
          <p className="mt-2 text-sm text-gray-500">
            Try: c1ccccc1 (benzene), CC(C)CC1=CC=C(C=C1)C(C)C(=O)O (ibuprofen)
          </p>
        </div>

        {pdbData && (
          <div className="bg-white p-6 rounded-lg shadow-md">
            <h2 className="text-2xl font-semibold mb-4 text-gray-800">
              3D Structure
            </h2>
            <div className="flex justify-center">
              <MoleculeViewer pdbData={pdbData} width={800} height={600} />
            </div>
            <p className="mt-4 text-sm text-gray-600 text-center">
              Use mouse to rotate • Scroll to zoom • Right-click to pan
            </p>
          </div>
        )}

        {!pdbData && !loading && (
          <div className="bg-white p-12 rounded-lg shadow-md text-center text-gray-500">
            <p className="text-lg">Enter a SMILES string above and click &quot;Visualize&quot; to see the 3D structure</p>
          </div>
        )}
      </div>
    </div>
  )
}
