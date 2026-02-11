import { useState, useEffect, useCallback } from "react"
import { ArrowLeft, FolderOpen } from "lucide-react"
import { cn } from "@/lib/utils"
import { Input } from "../ui/input"
import { AddWorkspaceContainer, AddWorkspaceStepHeader, AddWorkspacePrimaryButton } from "./primitives"
import { useT, use$T } from "@/context/LocaleContext"

interface AddWorkspaceStep_CreateNewProps {
  onBack: () => void
  onCreate: (folderPath: string, name: string) => Promise<void>
  isCreating: boolean
  /** 选中的应用名称，用于生成默认工作区名称 */
  appName?: string
}

/**
 * AddWorkspaceStep_CreateNew - 创建新工作区
 *
 * 流程：
 * 1. 必须选择文件夹
 * 2. 工作区名称默认为 "目录名-应用名"，选择目录后自动更新
 * 3. 创建时校验同名工作区
 */
export function AddWorkspaceStep_CreateNew({
  onBack,
  onCreate,
  isCreating,
  appName = ''
}: AddWorkspaceStep_CreateNewProps) {
  const t = useT()
  const $t = use$T()
  const [name, setName] = useState('')
  const [selectedPath, setSelectedPath] = useState<string | null>(null)
  const [error, setError] = useState<string | null>(null)
  const [isValidating, setIsValidating] = useState(false)
  // 标记用户是否手动修改过名称
  const [nameManuallyEdited, setNameManuallyEdited] = useState(false)

  // 根据选中的目录和应用名生成默认工作区名称
  const generateDefaultName = useCallback((folderPath: string, app: string) => {
    const folderName = folderPath.split(/[\\/]/).pop() || ''
    if (app && folderName) {
      return `${folderName}-${app}`
    }
    return folderName || app || ''
  }, [])

  // 校验工作区名称是否已存在
  useEffect(() => {
    const trimmedName = name.trim()
    if (!trimmedName) {
      setError(null)
      return
    }

    const validateName = async () => {
      setIsValidating(true)
      try {
        const result = await window.electronAPI.checkWorkspaceName(trimmedName)
        if (result.exists) {
          setError($t('名为 "{name}" 的工作区已存在，请修改名称', { name: trimmedName }))
        } else {
          setError(null)
        }
      } catch (err) {
        console.error('Failed to validate workspace name:', err)
      } finally {
        setIsValidating(false)
      }
    }

    // 防抖校验
    const timeout = setTimeout(validateName, 300)
    return () => clearTimeout(timeout)
  }, [name, $t])

  const handleBrowse = useCallback(async () => {
    const path = await window.electronAPI.openFolderDialog()
    if (path) {
      setSelectedPath(path)
      // 如果用户没有手动修改过名称，自动更新
      if (!nameManuallyEdited) {
        setName(generateDefaultName(path, appName))
      }
    }
  }, [appName, nameManuallyEdited, generateDefaultName])

  const handleNameChange = useCallback((e: React.ChangeEvent<HTMLInputElement>) => {
    setName(e.target.value)
    setNameManuallyEdited(true)
  }, [])

  const handleCreate = useCallback(async () => {
    if (!name.trim() || !selectedPath || error) return
    await onCreate(selectedPath, name.trim())
  }, [name, selectedPath, error, onCreate])

  const canCreate = name.trim() && selectedPath && !error && !isValidating && !isCreating

  return (
    <AddWorkspaceContainer>
      {/* 返回按钮 */}
      <button
        onClick={onBack}
        disabled={isCreating}
        className={cn(
          "self-start flex items-center gap-1 text-sm text-muted-foreground",
          "hover:text-foreground transition-colors mb-4",
          isCreating && "opacity-50 cursor-not-allowed"
        )}
      >
        <ArrowLeft className="h-4 w-4" />
        {t('返回')}
      </button>

      <AddWorkspaceStepHeader
        title={t('创建工作区')}
        description={t('选择工作区文件夹并输入名称。')}
      />

      <div className="mt-6 w-full space-y-6">
        {/* 选择文件夹 - 整行可点击 */}
        <div className="space-y-2">
          <label className="block text-sm font-medium text-foreground mb-2.5">
            {t('工作区文件夹')}
          </label>
          <button
            type="button"
            onClick={handleBrowse}
            disabled={isCreating}
            className={cn(
              "flex items-center gap-3 w-full p-3.5 rounded-xl text-left",
              "bg-background shadow-minimal",
              "border border-transparent",
              "transition-all duration-150",
              "hover:bg-foreground/5 hover:border-border/50",
              "focus:outline-none focus-visible:ring-2 focus-visible:ring-ring focus-visible:ring-offset-2",
              isCreating && "opacity-50 cursor-not-allowed"
            )}
          >
            <div className={cn(
              "flex h-9 w-9 shrink-0 items-center justify-center rounded-lg",
              selectedPath ? "bg-accent/10 text-accent" : "bg-foreground/5 text-muted-foreground"
            )}>
              <FolderOpen className="h-4.5 w-4.5" />
            </div>
            <div className="flex-1 min-w-0">
              {selectedPath ? (
                <p className="text-sm text-foreground truncate">{selectedPath}</p>
              ) : (
                <p className="text-sm text-muted-foreground">{t('点击选择文件夹')}</p>
              )}
            </div>
          </button>
        </div>

        {/* 工作区名称 */}
        <div className="space-y-2">
          <label className="block text-sm font-medium text-foreground mb-2.5">
            {t('工作区名称')}
          </label>
          <div className="bg-background shadow-minimal rounded-lg">
            <Input
              value={name}
              onChange={handleNameChange}
              placeholder={t('选择文件夹后自动生成')}
              disabled={isCreating}
              className="border-0 bg-transparent shadow-none"
            />
          </div>
          {error && (
            <p className="text-xs text-destructive">{error}</p>
          )}
        </div>

        {/* 创建按钮 */}
        <AddWorkspacePrimaryButton
          onClick={handleCreate}
          disabled={!canCreate}
          loading={isCreating}
          loadingText={t('创建中...')}
        >
          {t('创建')}
        </AddWorkspacePrimaryButton>
      </div>
    </AddWorkspaceContainer>
  )
}
