import craftLogo from "@/assets/craft_logo_c.svg"

interface CreatorFlowAppIconProps {
  className?: string
  size?: number
}

/**
 * CreatorFlowAppIcon - Displays the Craft logo (colorful "C" icon)
 */
export function CreatorFlowAppIcon({ className, size = 64 }: CreatorFlowAppIconProps) {
  return (
    <img
      src={craftLogo}
      alt="Craft"
      width={size}
      height={size}
      className={className}
    />
  )
}
