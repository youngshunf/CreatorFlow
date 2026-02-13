import { app } from './app'

const port = parseInt(process.env.PORT || '3000', 10)

console.log(`
╔═══════════════════════════════════════════════════════════╗
║                  Sprouty AI Cloud Server                    ║
╠═══════════════════════════════════════════════════════════╣
║  🚀 Server running on http://localhost:${port}               ║
║  📝 Health check: http://localhost:${port}/health            ║
║  🔧 Environment: ${process.env.NODE_ENV || 'development'}                         ║
╚═══════════════════════════════════════════════════════════╝
`)

export default {
  port,
  fetch: app.fetch,
}
