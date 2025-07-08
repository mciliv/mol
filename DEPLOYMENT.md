# Serverless Deployment Options

Your molecular analysis app is now configured for multiple serverless deployment options with minimal server management and dependencies.

## 🚀 Quick Deploy Options

### 1. **Vercel** (Recommended)
- ✅ Zero-config serverless functions
- ✅ Automatic HTTPS
- ✅ Global CDN
- ✅ Perfect for Node.js APIs

```bash
# Install Vercel CLI
npm install -g vercel

# Deploy
npm run deploy:vercel
```

### 2. **Netlify**
- ✅ Serverless functions
- ✅ Form handling
- ✅ Edge functions
- ✅ Great for static sites with APIs

```bash
# Install Netlify CLI
npm install -g netlify-cli

# Deploy
npm run deploy:netlify
```

### 3. **Railway** 
- ✅ Container-based deployment
- ✅ Database hosting
- ✅ Simple Git-based deployment

```bash
# Install Railway CLI
npm install -g @railway/cli

# Deploy
npm run deploy:railway
```

## 🔧 Environment Variables

Set these in your deployment platform:

```
OPENAI_API_KEY=your_openai_key
NODE_ENV=production
PORT=3000
```

## 🎯 Benefits of Serverless

- **No server management** - Platform handles scaling
- **Pay per use** - Only pay for actual usage
- **Auto-scaling** - Handles traffic spikes automatically
- **Zero downtime** - Automatic deployments
- **Global distribution** - CDN included

## 📱 Mobile-Friendly Features

Your app includes:
- HTTPS certificate generation
- Mobile camera access
- Cross-platform compatibility
- Responsive design

## 🛠️ One-Click Deploy

Click to deploy instantly:

[![Deploy to Vercel](https://vercel.com/button)](https://vercel.com/new/git/external?repository-url=https://github.com/mciliv/mol)
[![Deploy to Netlify](https://www.netlify.com/img/deploy/button.svg)](https://app.netlify.com/start/deploy?repository=https://github.com/mciliv/mol) 