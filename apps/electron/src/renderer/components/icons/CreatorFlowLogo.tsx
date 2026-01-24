import logoHorizontal from "@/assets/logo-horizontal.svg"

interface CreatorFlowLogoProps {
  className?: string
  height?: number
}

/**
 * CreatorFlow 横版 Logo - 图标 + 中英文名称
 */
export function CreatorFlowLogo({ className, height = 40 }: CreatorFlowLogoProps) {
  return (
    <img
      src={logoHorizontal}
      alt="创流 CreatorFlow"
      height={height}
      className={className}
    />
  )
}
