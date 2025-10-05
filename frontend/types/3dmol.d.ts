declare module '3dmol/build/3Dmol.js' {
  interface ViewerConfig {
    backgroundColor?: string
    defaultcolors?: Record<string, unknown>
  }

  interface ViewerGL {
    addModel(data: string, format: string): void
    setStyle(sel: Record<string, unknown>, style: Record<string, unknown>): void
    zoomTo(): void
    render(): void
    clear(): void
  }

  const $3Dmol: {
    createViewer(element: HTMLElement, config?: ViewerConfig): ViewerGL
  }

  export default $3Dmol
}
