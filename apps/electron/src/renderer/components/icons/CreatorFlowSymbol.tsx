interface CreatorFlowSymbolProps {
  className?: string
  size?: number
}

/**
 * 智小芽图标 - 双叶片种子设计 (青紫渐变)
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
        <linearGradient id="symbolGrad" x1="0%" y1="0%" x2="100%" y2="100%">
          <stop offset="0%" stopColor="#06B6D4"/>
          <stop offset="100%" stopColor="#8B5CF6"/>
        </linearGradient>
        <linearGradient id="symbolGrad2" x1="100%" y1="0%" x2="0%" y2="100%">
          <stop offset="0%" stopColor="#06B6D4"/>
          <stop offset="100%" stopColor="#8B5CF6"/>
        </linearGradient>
        <filter id="symbolGlow" x="-50%" y="-50%" width="200%" height="200%">
          <feGaussianBlur stdDeviation="1" result="blur"/>
          <feMerge>
            <feMergeNode in="blur"/>
            <feMergeNode in="SourceGraphic"/>
          </feMerge>
        </filter>
      </defs>

      {/* 圆形背景 */}
      <circle cx="32" cy="32" r="30" fill="#0F172A"/>

      {/* 内圈底座 */}
      <circle cx="32" cy="37" r="17" fill="#1E293B"/>

      {/* 左叶片 */}
      <path
        d="M30 36 Q22 30 20 17 Q25 25 29 32 Q30 35 30 36"
        fill="url(#symbolGrad)"
        filter="url(#symbolGlow)"
      />
      {/* 右叶片 */}
      <path
        d="M34 36 Q42 30 44 17 Q39 25 35 32 Q34 35 34 36"
        fill="url(#symbolGrad2)"
        filter="url(#symbolGlow)"
      />

      {/* 底部发光点 */}
      <circle cx="32" cy="43" r="3" fill="url(#symbolGrad)" filter="url(#symbolGlow)"/>
    </svg>
  )
}
