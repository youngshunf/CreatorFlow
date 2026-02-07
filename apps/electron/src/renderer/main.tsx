import React, { useState, useEffect } from 'react'
import ReactDOM from 'react-dom/client'
import { Provider as JotaiProvider } from 'jotai'
import App from './App'
import { ThemeProvider } from './context/ThemeContext'
import { LocaleProvider } from './context/LocaleContext'
import { Toaster } from '@/components/ui/sonner'
import './index.css'

/**
 * Root component - loads workspace ID for theme context and renders App
 * App.tsx handles window mode detection internally (main vs tab-content)
 */
function Root() {
  // Load workspace ID for theme context (workspace-specific theme overrides)
  const [workspaceId, setWorkspaceId] = useState<string | null>(null)

  useEffect(() => {
    window.electronAPI?.getWindowWorkspace?.().then((id) => {
      setWorkspaceId(id)
    })
  }, [])

  return (
    <ThemeProvider activeWorkspaceId={workspaceId}>
      <App />
      <Toaster />
    </ThemeProvider>
  )
}

ReactDOM.createRoot(document.getElementById('root')!).render(
  <React.StrictMode>
    <JotaiProvider>
      <ThemeProvider>
        <LocaleProvider>
          <Root />
          <Toaster />
        </LocaleProvider>
      </ThemeProvider>
    </JotaiProvider>
  </React.StrictMode>
)
