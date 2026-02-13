import logoHorizontal from "@/assets/logo-horizontal.svg"

interface SproutyLogoProps {
  className?: string
  height?: number
}

/**
 * 智小芽横版 Logo - 图标 + 中文名称
 */
export function SproutyLogo({ className, height = 40 }: SproutyLogoProps) {
  return (
    <img
      src={logoHorizontal}
      alt="智小芽"
      height={height}
      className={className}
    />
  )
}
