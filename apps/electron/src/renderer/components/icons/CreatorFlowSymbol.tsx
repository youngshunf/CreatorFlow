interface CreatorFlowSymbolProps {
  className?: string
  size?: number
}

/**
 * CreatorFlow 图标 - 流动曲线 + 笔触设计
 */
export function CreatorFlowSymbol({ className, size = 32 }: CreatorFlowSymbolProps) {
  return (
    <svg
      width={size}
      height={size}
      viewBox="0 0 64 64"
      className={className}
      fill="none"
      xmlns="http://www.w3.org/2000/svg"
    >
      <defs>
        <linearGradient id="iconGradient" x1="0%" y1="0%" x2="100%" y2="100%">
          <stop offset="0%" stopColor="#3B82F6"/>
          <stop offset="50%" stopColor="#6366F1"/>
          <stop offset="100%" stopColor="#8B5CF6"/>
        </linearGradient>
      </defs>

      {/* 圆角矩形背景 */}
      <rect x="5" y="5" width="54" height="54" rx="13" fill="url(#iconGradient)"/>

      {/* 主流动曲线 */}
      <path
        d="M16 45 C21 38, 24 32, 27 30 S34 28, 37 27 S45 22, 48 19"
        stroke="white"
        strokeWidth="3.5"
        strokeLinecap="round"
        fill="none"
        opacity="0.9"
      />

      {/* 辅助流动线 */}
      <path
        d="M19 48 C24 42, 27 37, 30 35 S37 32, 40 30"
        stroke="white"
        strokeWidth="1.5"
        strokeLinecap="round"
        fill="none"
        opacity="0.4"
      />

      {/* 创作笔触 */}
      <g transform="translate(18, 16)">
        <path d="M0 11 L5 6 L6.5 7.5 L1.5 12.5 Z" fill="white" opacity="0.95"/>
        <path d="M5 6 L11 0" stroke="white" strokeWidth="2" strokeLinecap="round" opacity="0.9"/>
        <path d="M1.5 12.5 L0 16" stroke="white" strokeWidth="1.8" strokeLinecap="round" opacity="0.8"/>
      </g>

      {/* 灵感粒子 */}
      <circle cx="48" cy="19" r="3.2" fill="white" opacity="0.95"/>
      <circle cx="53" cy="25" r="1.6" fill="white" opacity="0.6"/>
      <circle cx="45" cy="14" r="1.3" fill="white" opacity="0.45"/>
      <circle cx="50" cy="30" r="1" fill="white" opacity="0.45"/>
    </svg>
  )
}
