import logoIcon from "@/assets/logo-icon.svg"

interface CreatorFlowAppIconProps {
  className?: string
  size?: number
}

/**
 * CreatorFlowAppIcon - 显示智小芽 App 图标
 */
export function CreatorFlowAppIcon({ className, size = 64 }: CreatorFlowAppIconProps) {
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
