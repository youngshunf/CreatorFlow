import React, { useState, useEffect } from 'react'
import ReactDOM from 'react-dom/client'
import { Provider as JotaiProvider } from 'jotai'
import App from './App'
import { ThemeProvider } from './context/ThemeContext'
import { LocaleProvider } from './context/LocaleContext'
import { Toaster } from '@/components/ui/sonner'
import './index.css'

/**
 * 应用崩溃时的最小回退 UI
 */
function CrashFallback() {
  return (
    <div className="flex flex-col items-center justify-center h-screen font-sans text-foreground/50 gap-3">
      <p className="text-base font-medium">Something went wrong</p>
      <p className="text-[13px]">Please restart the app. The error has been reported.</p>
      <button
        onClick={() => window.location.reload()}
        className="mt-2 px-4 py-1.5 rounded-md bg-background shadow-minimal text-[13px] text-foreground/70 cursor-pointer"
      >
        Reload
      </button>
    </div>
  )
}

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
