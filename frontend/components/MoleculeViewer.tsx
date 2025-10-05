'use client'

import dynamic from 'next/dynamic'
import { useEffect, useRef, useCallback } from 'react'

interface MoleculeViewerProps {
  pdbData: string
  width?: number
  height?: number
}

const MoleculeViewer: React.FC<MoleculeViewerProps> = ({
  pdbData,
  width = 800,
  height = 600
}) => {
  const viewerRef = useRef<HTMLDivElement>(null)
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const $3Dmol = useRef<any>(null)
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const viewerInstance = useRef<any>(null)

  const initializeViewer = useCallback(() => {
    if (!viewerRef.current || !$3Dmol.current) return

    // Clear existing viewer
    if (viewerInstance.current) {
      viewerInstance.current.clear()
    }

    const viewer = $3Dmol.current.createViewer(viewerRef.current, {
      backgroundColor: 'white',
    })

    viewer.addModel(pdbData, 'pdb')
    viewer.setStyle({}, { stick: {} })
    viewer.zoomTo()
    viewer.render()

    viewerInstance.current = viewer
  }, [pdbData])

  useEffect(() => {
    // Load 3Dmol.js only on client
    if (typeof window !== 'undefined' && !$3Dmol.current) {
      import('3dmol/build/3Dmol.js').then((module) => {
        $3Dmol.current = module.default
        initializeViewer()
      })
    }
  }, [initializeViewer])

  useEffect(() => {
    // Update viewer when pdbData changes
    if ($3Dmol.current && pdbData) {
      initializeViewer()
    }
  }, [pdbData, initializeViewer])

  return (
    <div
      ref={viewerRef}
      style={{ width: `${width}px`, height: `${height}px`, position: 'relative' }}
      className="border border-gray-300 rounded-lg shadow-sm"
    />
  )
}

// Export with SSR disabled for WebGL
export default dynamic(() => Promise.resolve(MoleculeViewer), {
  ssr: false, // Critical: disable SSR for WebGL
})
