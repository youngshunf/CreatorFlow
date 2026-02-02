import logoHorizontal from "@/assets/logo-horizontal.svg"

interface CreatorFlowLogoProps {
  className?: string
  height?: number
}

/**
 * 智小芽横版 Logo - 图标 + 中文名称
 */
export function CreatorFlowLogo({ className, height = 40 }: CreatorFlowLogoProps) {
  return (
    <img
      src={logoHorizontal}
      alt="智小芽"
      height={height}
      className={className}
    />
  )
}
