import { useState, useEffect, useCallback } from "react"
import { ArrowLeft } from "lucide-react"
import { cn } from "@/lib/utils"
import { slugify } from "@/lib/slugify"
import { Input } from "../ui/input"
import { Button } from "../ui/button"
import { AddWorkspaceContainer, AddWorkspaceStepHeader, AddWorkspaceSecondaryButton, AddWorkspacePrimaryButton } from "./primitives"
import { AddWorkspace_RadioOption } from "./AddWorkspace_RadioOption"
import { useT, use$T } from "@/context/LocaleContext"

type LocationOption = 'default' | 'custom'

interface AddWorkspaceStep_CreateNewProps {
  onBack: () => void
  onCreate: (folderPath: string, name: string) => Promise<void>
  isCreating: boolean
}

/**
 * AddWorkspaceStep_CreateNew - Create a new workspace
 *
 * Fields:
 * - Workspace name (required)
 * - Location: Default (~/.creator-flow/workspaces/) or Custom
 */
export function AddWorkspaceStep_CreateNew({
  onBack,
  onCreate,
  isCreating
}: AddWorkspaceStep_CreateNewProps) {
  const t = useT()
  const $t = use$T()
  const [name, setName] = useState('')
  const [locationOption, setLocationOption] = useState<LocationOption>('default')
  const [customPath, setCustomPath] = useState<string | null>(null)
  const [homeDir, setHomeDir] = useState('')
  const [error, setError] = useState<string | null>(null)
  const [isValidating, setIsValidating] = useState(false)

  // Get home directory on mount
  useEffect(() => {
    window.electronAPI.getHomeDir().then(setHomeDir)
  }, [])

  const slug = slugify(name)
  const defaultBasePath = homeDir ? `${homeDir}/.creator-flow/workspaces` : '~/.creator-flow/workspaces'
  const finalPath = locationOption === 'default'
    ? `${defaultBasePath}/${slug}`
    : customPath
      ? `${customPath}/${slug}`
      : null

  // Validate slug uniqueness when name changes
  useEffect(() => {
    if (!slug) {
      setError(null)
      return
    }

    const validateSlug = async () => {
      setIsValidating(true)
      try {
        const result = await window.electronAPI.checkWorkspaceSlug(slug)
        if (result.exists) {
          setError($t('名为 "{name}" 的工作区已存在', { name: slug }))
        } else {
          setError(null)
        }
      } catch (err) {
        console.error('Failed to validate workspace slug:', err)
      } finally {
        setIsValidating(false)
      }
    }

    // Debounce validation
    const timeout = setTimeout(validateSlug, 300)
    return () => clearTimeout(timeout)
  }, [slug])

  const handleBrowse = useCallback(async () => {
    const path = await window.electronAPI.openFolderDialog()
    if (path) {
      setCustomPath(path)
    }
  }, [])

  const handleCreate = useCallback(async () => {
    if (!name.trim() || !finalPath || error) return
    await onCreate(finalPath, name.trim())
  }, [name, finalPath, error, onCreate])

  const canCreate = name.trim() && finalPath && !error && !isValidating && !isCreating

  return (
    <AddWorkspaceContainer>
      {/* Back button */}
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
        description={t('输入名称并选择工作区的存储位置。')}
      />

      <div className="mt-6 w-full space-y-6">
        {/* Workspace name */}
        <div className="space-y-2">
          <label className="block text-sm font-medium text-foreground mb-2.5">
            {t('工作区名称')}
          </label>
          <div className="bg-background shadow-minimal rounded-lg">
            <Input
              value={name}
              onChange={(e) => setName(e.target.value)}
              placeholder={t('我的工作区')}
              disabled={isCreating}
              autoFocus
              className="border-0 bg-transparent shadow-none"
            />
          </div>
          {error && (
            <p className="text-xs text-destructive">{error}</p>
          )}
        </div>

        {/* Location selection */}
        <div className="space-y-3">
          <label className="block text-sm font-medium text-foreground mb-2.5">
            {t('位置')}
          </label>

          {/* Default location option */}
          <AddWorkspace_RadioOption
            name="location"
            checked={locationOption === 'default'}
            onChange={() => setLocationOption('default')}
            disabled={isCreating}
            title={t('默认位置')}
            subtitle={t('在 .creator-flow 文件夹下')}
          />

          {/* Custom location option */}
          <AddWorkspace_RadioOption
            name="location"
            checked={locationOption === 'custom'}
            onChange={() => setLocationOption('custom')}
            disabled={isCreating}
            title={t('选择位置')}
            subtitle={customPath || t('选择新工作区的存储位置。')}
            action={locationOption === 'custom' ? (
              <AddWorkspaceSecondaryButton
                onClick={(e) => {
                  e.preventDefault()
                  handleBrowse()
                }}
                disabled={isCreating}
              >
                {t('浏览')}
              </AddWorkspaceSecondaryButton>
            ) : undefined}
          />
        </div>

        {/* Create button */}
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
