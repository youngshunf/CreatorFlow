import { useState, useEffect, useCallback, useMemo } from "react"
import { X, AlertTriangle } from "lucide-react"
import { motion } from "motion/react"
import { Dithering } from "@paper-design/shaders-react"
import { FullscreenOverlayBase } from "@sprouty-ai/ui"
import { cn } from "@/lib/utils"
import { overlayTransitionIn } from "@/lib/animations"
import { AddWorkspaceStep_ChooseApp } from "./AddWorkspaceStep_ChooseApp"
import { AddWorkspaceStep_CreateNew } from "./AddWorkspaceStep_CreateNew"
import { Button } from "@/components/ui/button"
import { useT } from "@/context/LocaleContext"
import type { Workspace } from "../../../shared/types"

type CreationStep = 'choose-app' | 'create'
type AppSource = 'bundled' | 'marketplace'

interface WorkspaceCreationScreenProps {
  /** Callback when a workspace is created successfully */
  onWorkspaceCreated: (workspace: Workspace) => void
  /** Callback when the screen is dismissed */
  onClose: () => void
  /** Initial marketplace app to use (skips app selection step) */
  initialMarketplaceApp?: { id: string; name: string }
  className?: string
}

/**
 * WorkspaceCreationScreen - Full-screen overlay for creating workspaces
 *
 * Flow:
 * 1. Choose App: Select application template
 * 2. Create: Enter name + choose location
 */
export function WorkspaceCreationScreen({
  onWorkspaceCreated,
  onClose,
  initialMarketplaceApp,
  className
}: WorkspaceCreationScreenProps) {
  const t = useT()
  // If initialMarketplaceApp is provided, skip to create step
  const [step, setStep] = useState<CreationStep>(
    initialMarketplaceApp ? 'create' : 'choose-app'
  )
  const [selectedAppId, setSelectedAppId] = useState<string>(
    initialMarketplaceApp ? `marketplace:${initialMarketplaceApp.id}` : 'app.general'
  )
  const [selectedAppName, setSelectedAppName] = useState<string>(
    initialMarketplaceApp?.name || ''
  )
  const [selectedAppSource, setSelectedAppSource] = useState<AppSource>(
    initialMarketplaceApp ? 'marketplace' : 'bundled'
  )
  const [isCreating, setIsCreating] = useState(false)
  const [creatingStage, setCreatingStage] = useState<string>('')
  const [dimensions, setDimensions] = useState({ width: 1920, height: 1080 })
  
  // State for existing app confirmation dialog
  const [showExistingAppDialog, setShowExistingAppDialog] = useState(false)
  const [existingAppInfo, setExistingAppInfo] = useState<{ name: string; version: string } | null>(null)
  const [pendingCreateParams, setPendingCreateParams] = useState<{ folderPath: string; name: string } | null>(null)

  // Track window dimensions for shader
  useEffect(() => {
    const updateDimensions = () => {
      setDimensions({ width: window.innerWidth, height: window.innerHeight })
    }
    updateDimensions()
    window.addEventListener('resize', updateDimensions)
    return () => window.removeEventListener('resize', updateDimensions)
  }, [])

  // Wrap onClose to prevent closing during creation
  // FullscreenOverlayBase handles ESC key, this wrapper prevents closing when busy
  const handleClose = useCallback(() => {
    if (!isCreating) {
      onClose()
    }
  }, [isCreating, onClose])

  const handleChooseApp = useCallback((appId: string, appName: string, source: AppSource) => {
    setSelectedAppId(appId)
    setSelectedAppName(appName)
    setSelectedAppSource(source)
    setStep('create')
  }, [])

  const handleCreateWorkspace = useCallback(async (folderPath: string, name: string, installMode?: 'force' | 'merge') => {
    setIsCreating(true)
    setCreatingStage(t('正在创建工作区...'))
    try {
      // For marketplace apps, extract the actual app ID
      const appId = selectedAppSource === 'marketplace'
        ? selectedAppId.replace('marketplace:', '')
        : selectedAppId

      // Show different stages for marketplace apps
      if (selectedAppSource === 'marketplace') {
        setCreatingStage(t('正在下载应用...'))
        // Simulate stage updates (since we don't have real progress events yet)
        setTimeout(() => setCreatingStage(t('正在安装技能...')), 2000)
        setTimeout(() => setCreatingStage(t('正在初始化工作区...')), 4000)
      }

      // Pass the selected app ID and source to workspace creation
      const result = await window.electronAPI.createWorkspace(
        folderPath,
        name,
        appId,
        selectedAppSource === 'marketplace' ? 'marketplace' : undefined,
        installMode
      )

      console.log('[WorkspaceCreation] createWorkspace result:', result)

      // Check if there's an existing app conflict (workspace will be null)
      if (result.existingApp && !result.workspace) {
        console.log('[WorkspaceCreation] Showing existing app dialog:', result.existingApp)
        // Show confirmation dialog
        setExistingAppInfo(result.existingApp)
        setPendingCreateParams({ folderPath, name })
        setShowExistingAppDialog(true)
        return
      }

      if (result.workspace) {
        setCreatingStage(t('完成！'))
        onWorkspaceCreated(result.workspace)
      }
    } finally {
      setIsCreating(false)
      setCreatingStage('')
    }
  }, [onWorkspaceCreated, selectedAppId, selectedAppSource, t])

  // Handle force install (overwrite with backup) after user confirms
  const handleForceInstall = useCallback(async () => {
    if (!pendingCreateParams) return
    setShowExistingAppDialog(false)
    await handleCreateWorkspace(pendingCreateParams.folderPath, pendingCreateParams.name, 'force')
  }, [pendingCreateParams, handleCreateWorkspace])

  // Handle merge install (keep existing, add new files)
  const handleMergeInstall = useCallback(async () => {
    if (!pendingCreateParams) return
    setShowExistingAppDialog(false)
    await handleCreateWorkspace(pendingCreateParams.folderPath, pendingCreateParams.name, 'merge')
  }, [pendingCreateParams, handleCreateWorkspace])

  // Handle choosing different folder
  const handleChooseDifferentFolder = useCallback(() => {
    setShowExistingAppDialog(false)
    setExistingAppInfo(null)
    setPendingCreateParams(null)
    // Stay on create step so user can choose different folder
  }, [])

  const renderStep = () => {
    switch (step) {
      case 'choose-app':
        return (
          <AddWorkspaceStep_ChooseApp
            onNext={handleChooseApp}
            isLoading={isCreating}
            isFirstStep={true}
          />
        )

      case 'create':
        return (
          <AddWorkspaceStep_CreateNew
            onBack={() => setStep('choose-app')}
            onCreate={handleCreateWorkspace}
            isCreating={isCreating}
            defaultName={selectedAppName}
          />
        )

      default:
        return null
    }
  }

  // Get theme colors from CSS variables for the shader
  const shaderColors = useMemo(() => {
    if (typeof window === 'undefined') return { back: '#00000000', front: '#684e85' }
    const root = document.documentElement
    const isDark = root.classList.contains('dark')
    // Transparent back, accent-tinted front
    return isDark
      ? { back: '#00000000', front: '#9b7bb8' }  // lighter accent for dark mode
      : { back: '#00000000', front: '#684e85' }  // accent color
  }, [])

  // FullscreenOverlayBase handles portal, traffic lights, and ESC key
  return (
    <FullscreenOverlayBase
      isOpen={true}
      onClose={handleClose}
      className={cn("z-splash flex flex-col bg-background", className)}
    >
      <motion.div
        initial={{ opacity: 0 }}
        animate={{ opacity: 1 }}
        exit={{ opacity: 0 }}
        transition={overlayTransitionIn}
        className="flex flex-col flex-1"
      >
        {/* Dithering shader background */}
        <motion.div
          initial={{ opacity: 0 }}
          animate={{ opacity: 0.3 }}
          transition={overlayTransitionIn}
          className="absolute inset-0 pointer-events-none"
        >
          <Dithering
            colorBack={shaderColors.back}
            colorFront={shaderColors.front}
            shape="swirl"
            type="8x8"
            size={2}
            speed={1}
            scale={1}
            width={dimensions.width}
            height={dimensions.height}
          />
        </motion.div>

        {/* Header with drag region and close button */}
        <header className="titlebar-drag-region relative h-[50px] shrink-0 flex items-center justify-end px-6">
          {/* Close button - explicitly no-drag */}
          <motion.button
            initial={{ opacity: 0 }}
            animate={{ opacity: 1 }}
            transition={overlayTransitionIn}
            onClick={(e) => {
              e.stopPropagation()
              handleClose()
            }}
            disabled={isCreating}
            className={cn(
              "titlebar-no-drag flex items-center justify-center p-2 rounded-[6px]",
              "bg-background shadow-minimal hover:bg-foreground-5",
              "text-muted-foreground hover:text-foreground",
              "transition-colors focus:outline-none focus-visible:ring-2 focus-visible:ring-ring",
              "mr-[-8px] mt-2",
              isCreating && "opacity-50 cursor-not-allowed"
            )}
            aria-label="Close"
          >
            <X className="h-4 w-4" />
          </motion.button>
        </header>

        {/* Main content */}
        <motion.main
          initial={{ opacity: 0 }}
          animate={{ opacity: 1 }}
          transition={overlayTransitionIn}
          className={cn(
            "relative flex flex-1 overflow-hidden",
            step === 'choose-app' ? 'p-0' : 'items-center justify-center p-8'
          )}
        >
          {renderStep()}
        </motion.main>

        {/* Existing App Confirmation Dialog */}
        {showExistingAppDialog && existingAppInfo && (
          <div className="absolute inset-0 z-50 flex items-center justify-center bg-background/80 backdrop-blur-sm">
            <motion.div
              initial={{ opacity: 0, scale: 0.95 }}
              animate={{ opacity: 1, scale: 1 }}
              className="w-full max-w-md rounded-xl border bg-background p-6 shadow-lg"
            >
              <div className="flex items-start gap-4">
                <div className="flex h-10 w-10 shrink-0 items-center justify-center rounded-full bg-yellow-500/10">
                  <AlertTriangle className="h-5 w-5 text-yellow-500" />
                </div>
                <div className="flex-1">
                  <h3 className="text-lg font-semibold">
                    {t('工作区已有应用')}
                  </h3>
                  <p className="mt-2 text-sm text-muted-foreground">
                    {t('该工作区已安装应用')} <span className="font-medium text-foreground">"{existingAppInfo.name}"</span> ({existingAppInfo.version})。
                  </p>
                  <p className="mt-2 text-sm text-muted-foreground">
                    {t('您可以选择其他工作区，或选择安装方式：')}
                  </p>
                  <ul className="mt-2 text-sm text-muted-foreground list-disc list-inside space-y-1">
                    <li>{t('合并安装：保留现有文件，只添加/更新新文件')}</li>
                    <li>{t('覆盖安装：先备份现有数据，再安装新应用')}</li>
                  </ul>
                </div>
              </div>
              <div className="mt-6 flex flex-col gap-2">
                <div className="flex gap-2">
                  <Button
                    variant="default"
                    className="flex-1"
                    onClick={handleMergeInstall}
                  >
                    {t('合并安装')}
                  </Button>
                  <Button
                    variant="outline"
                    className="flex-1"
                    onClick={handleForceInstall}
                  >
                    {t('覆盖安装')}
                  </Button>
                </div>
                <Button
                  variant="ghost"
                  className="w-full"
                  onClick={handleChooseDifferentFolder}
                >
                  {t('选择其他工作区')}
                </Button>
              </div>
            </motion.div>
          </div>
        )}

        {/* Creating Progress Overlay */}
        {isCreating && creatingStage && (
          <div className="absolute inset-0 z-40 flex items-center justify-center bg-background/80 backdrop-blur-sm">
            <motion.div
              initial={{ opacity: 0, scale: 0.95 }}
              animate={{ opacity: 1, scale: 1 }}
              className="w-full max-w-md rounded-xl border bg-background p-8 shadow-lg"
            >
              <div className="flex flex-col items-center gap-4">
                <div className="h-12 w-12 animate-spin rounded-full border-4 border-accent/20 border-t-accent" />
                <div className="text-center">
                  <h3 className="text-lg font-semibold">{creatingStage}</h3>
                  <p className="mt-2 text-sm text-muted-foreground">
                    {t('请稍候，这可能需要几分钟...')}
                  </p>
                </div>
              </div>
            </motion.div>
          </div>
        )}
      </motion.div>
    </FullscreenOverlayBase>
  )
}
