import logoIcon from "@/assets/logo-icon.svg"

interface SproutyAppIconProps {
  className?: string
  size?: number
}

/**
 * SproutyAppIcon - 显示智小芽 App 图标
 */
export function SproutyAppIcon({ className, size = 64 }: SproutyAppIconProps) {
  return (
    <img
      src={logoIcon}
      alt="智小芽"
      width={size}
      height={size}
      className={className}
    />
  )
}
