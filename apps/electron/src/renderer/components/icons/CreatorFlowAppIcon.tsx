import logoIcon from "@/assets/logo-icon.svg"

interface CreatorFlowAppIconProps {
  className?: string
  size?: number
}

/**
 * CreatorFlowAppIcon - 显示创流 App 图标
 */
export function CreatorFlowAppIcon({ className, size = 64 }: CreatorFlowAppIconProps) {
  return (
    <img
      src={logoIcon}
      alt="创流"
      width={size}
      height={size}
      className={className}
    />
  )
}
