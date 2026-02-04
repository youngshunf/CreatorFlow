import { motion } from 'motion/react'
import { CreatorFlowSymbol } from './icons/CreatorFlowSymbol'

interface SplashScreenProps {
  isExiting: boolean
  onExitComplete?: () => void
}

/**
 * SplashScreen - 智小芽启动页面
 *
 * 显示品牌信息、slogan 和加载进度
 */
export function SplashScreen({ isExiting, onExitComplete }: SplashScreenProps) {
  return (
    <motion.div
      className="fixed inset-0 z-splash flex flex-col items-center justify-center bg-background overflow-hidden"
      initial={{ opacity: 1 }}
      animate={{ opacity: isExiting ? 0 : 1 }}
      transition={{ duration: 0.5, ease: 'easeOut' }}
      onAnimationComplete={() => {
        if (isExiting && onExitComplete) {
          onExitComplete()
        }
      }}
    >
      {/* 品牌图标 */}
      <motion.div
        initial={{ scale: 0.8, opacity: 0 }}
        animate={{
          scale: isExiting ? 1.2 : 1,
          opacity: isExiting ? 0 : 1
        }}
        transition={{
          duration: isExiting ? 0.3 : 0.6,
          ease: [0.16, 1, 0.3, 1]
        }}
        className="mb-8"
      >
        <CreatorFlowSymbol className="h-16 text-accent" />
      </motion.div>

      {/* 品牌名称 */}
      <motion.h1
        initial={{ y: 20, opacity: 0 }}
        animate={{
          y: isExiting ? -10 : 0,
          opacity: isExiting ? 0 : 1
        }}
        transition={{
          duration: isExiting ? 0.3 : 0.6,
          delay: isExiting ? 0 : 0.2,
          ease: [0.16, 1, 0.3, 1]
        }}
        className="text-3xl font-bold text-foreground mb-3"
      >
        智小芽
      </motion.h1>

      {/* Slogan */}
      <motion.p
        initial={{ y: 20, opacity: 0 }}
        animate={{
          y: isExiting ? -10 : 0,
          opacity: isExiting ? 0 : 1
        }}
        transition={{
          duration: isExiting ? 0.3 : 0.6,
          delay: isExiting ? 0 : 0.3,
          ease: [0.16, 1, 0.3, 1]
        }}
        className="text-sm text-muted-foreground mb-16"
      >
        超级桌面AGENT，越用越懂你
      </motion.p>

      {/* 加载进度条 */}
      <motion.div
        initial={{ opacity: 0 }}
        animate={{ opacity: isExiting ? 0 : 1 }}
        transition={{
          duration: isExiting ? 0.2 : 0.4,
          delay: isExiting ? 0 : 0.4
        }}
        className="absolute bottom-12 left-1/2 -translate-x-1/2 w-64"
      >
        <div className="h-1 bg-foreground/10 rounded-full overflow-hidden">
          <motion.div
            className="h-full bg-accent rounded-full"
            initial={{ width: '0%' }}
            animate={{ width: isExiting ? '100%' : '60%' }}
            transition={{
              duration: isExiting ? 0.3 : 2,
              ease: 'easeInOut'
            }}
          />
        </div>
      </motion.div>
    </motion.div>
  )
}
