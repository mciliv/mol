# 🚀 How to Run Flutter Molecular Analyzer

## The Flutter Equivalent of "Local Development Server"

In Flutter, the equivalent of your web app's **local development server** is called:

### **`flutter run`** - Hot Reload Development Server

This command:
- ✅ Starts the app in development mode
- ✅ Enables **hot reload** (instant UI updates when you save files)
- ✅ Provides debug console and error reporting
- ✅ Automatically rebuilds when code changes

## 📦 Prerequisites

### 1. Install Flutter
```bash
# macOS (using Homebrew)
brew install --cask flutter

# Or download from: https://flutter.dev/docs/get-started/install
```

### 2. Verify Installation
```bash
flutter doctor
# This checks your Flutter installation and dependencies
```

## ▶️ Running the App

### 1. Navigate to Flutter Project
```bash
cd flutter_prototype
```

### 2. Get Dependencies
```bash
flutter pub get
# This installs all dependencies from pubspec.yaml
# Equivalent to "npm install" in web development
```

### 3. Run Development Server
```bash
# For Web (like your current setup)
flutter run -d web

# For Chrome specifically  
flutter run -d chrome --web-port 3000

# For mobile (if you have simulator/device)
flutter run

# For desktop
flutter run -d macos    # or windows/linux
```

## 🔥 Hot Reload Commands

While the app is running:
- **`r`** - Hot reload (updates UI instantly)
- **`R`** - Hot restart (full app restart)
- **`q`** - Quit

## 📱 Platform-Specific Commands

```bash
# Web development (equivalent to your current workflow)
flutter run -d web --web-port 3000 --web-hostname 0.0.0.0

# Mobile testing
flutter run -d ios          # iOS simulator
flutter run -d android      # Android emulator

# Desktop development
flutter run -d macos        # Native macOS app
flutter run -d windows      # Native Windows app
flutter run -d linux        # Native Linux app
```

## 🛠️ Development Workflow

### 1. Start Development Server
```bash
flutter run -d web
```

### 2. Edit Code
- Open any `.dart` file
- Make changes
- Save file
- **Instant hot reload!** ⚡

### 3. Debug
```bash
# Enable debug mode with DevTools
flutter run -d web --debug

# View logs
flutter logs
```

## 🔧 Build for Production

```bash
# Web build (equivalent to your production build)
flutter build web

# Mobile builds
flutter build apk          # Android
flutter build ipa          # iOS
flutter build macos        # macOS app
```

## 📂 Project Structure

```
flutter_prototype/
├── lib/
│   ├── main.dart           # App entry point
│   ├── screens/            # Page layouts  
│   ├── widgets/            # UI components
│   └── services/           # Business logic
├── pubspec.yaml            # Dependencies (like package.json)
└── web/                    # Web-specific files (auto-generated)
```

## 💡 Key Differences from Web Development

| **Web Development** | **Flutter Development** |
|-------------------|----------------------|
| `npm start` | `flutter run` |
| `npm install` | `flutter pub get` |
| `npm run build` | `flutter build web` |
| Browser DevTools | Flutter DevTools |
| HTML/CSS/JS | Single Dart language |
| Manual responsive design | Built-in adaptive widgets |

## 🎯 What You'll See

1. **Payment section** appears at top when no payment method
2. **Text input** for molecule names  
3. **Camera/Photo tabs** for image analysis
4. **Results area** showing molecular analysis
5. **Hot reload** - change any code and see instant updates!

## 🚧 Current Status

This is a **working prototype** that demonstrates:
- ✅ Full app structure and navigation
- ✅ Payment UI without modal complexity
- ✅ State management with Provider
- ✅ Camera integration setup
- ✅ Results display with 3D placeholders
- 🔄 Ready for real Stripe/API integration

**Ready to run once Flutter is installed!** 