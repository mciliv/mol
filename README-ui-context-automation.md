# UI Context Automation - Usage Guide

## 🎯 Automatic Context for Iterative UI Development

The UI Context Automation system automatically captures your UI state and generates formatted context for follow-up prompts, making iterative development more accurate and efficient.

## ✨ Features

- **Automatic UI State Capture**: Viewport, layout, CSS variables, responsive state
- **Rule Compliance Checking**: Validates against all .mdc rules (layout-system, product-requirements, etc.)
- **Overlap Detection**: Automatically detects and measures element overlaps
- **Change Tracking**: Records UI changes with timestamps
- **Formatted Context**: Ready-to-paste context for AI prompts
- **Development Integration**: Auto-enables on localhost

## 🚀 Quick Start

### Auto-Enabled (Localhost)
The system automatically activates when you load the app on localhost.

### Manual Control
```javascript
// Enable/disable
window.app.enableUIContextAutomation()

// Check current rule compliance
window.app.checkRuleCompliance()
// Returns: { compliant: true/false, issues: [...], score: 85 }

// Capture UI state snapshot
window.app.captureUIState('before_changes')

// Export formatted context for next prompt
window.app.exportContextForPrompt()
// Copies formatted context to clipboard
```

## 📋 Generated Context Format

When you run `window.app.exportContextForPrompt()`, you get:

```markdown
## Current UI State (Auto-Generated 2025-01-18T...)
**Viewport**: 1920x1080 (desktop)
**Rule Compliance**: 95%

**Active Issues** (2):
- [HIGH] compliance: Text input not using safe area margin
- [MEDIUM] overlap: account-link ↔ object-input (severity: 15.2%)

**Layout State**:
- Sidebar visible: false
- Safe area right: calc(300px + 20px)
- Focused element: object-input
- Visible sections: camera-container, results-section

**Overlapping Elements**: DETECTED
  - account-link ↔ object-input (severity: 15.2%)

**Suggested Focus**: Address element overlaps using layout debug mode

**Recent Changes** (3):
- 2025-01-18T12:30:45Z: DOM mutation detected
- 2025-01-18T12:30:40Z: Viewport resized to 1920x1080
- 2025-01-18T12:30:35Z: Card button position updated
```

## 🔍 What Gets Captured

### UI State
- **Viewport**: Width, height, device pixel ratio
- **Layout**: Sidebar visibility, element positions, safe areas
- **CSS Variables**: All custom properties (--space-*, --z-*, etc.)
- **Responsive State**: Current breakpoint (mobile/tablet/desktop)
- **Active Elements**: Focused element, visible sections, error states

### Rule Compliance
- **Layout System**: CSS custom properties usage, safe areas, z-index stack
- **Product Requirements**: Card button presence, borderless design, black backgrounds
- **UI Memory**: Prevention of known issues

### Overlap Detection
- **Element Analysis**: Checks key elements for overlaps
- **Severity Calculation**: Measures overlap percentage
- **Priority Rating**: High/medium/low based on severity

## 📝 Integration with AI Prompts

### Before Making Changes
```javascript
// Capture current state
window.app.captureUIState('before_fix')

// Check what needs fixing
window.app.checkRuleCompliance()
```

### After Making Changes  
```javascript
// Capture new state
window.app.captureUIState('after_fix')

// Generate context for next iteration
window.app.exportContextForPrompt()
// Paste this into your next AI prompt
```

### For Complex Issues
1. Enable layout debug: `window.app.enableLayoutDebug()`
2. Export context: `window.app.exportContextForPrompt()`
3. Include both in your AI prompt for comprehensive troubleshooting

## 🛠️ Development Workflow

### Recommended Pattern
1. **Start**: System auto-captures baseline state
2. **Plan**: Check rule compliance to identify issues
3. **Change**: Make your UI modifications
4. **Validate**: Export context to see what changed
5. **Iterate**: Use context in next AI prompt for refinement

### Debug Mode Integration
```javascript
// Visual debugging
window.app.enableLayoutDebug()  // Shows element boundaries

// Context debugging  
window.app.exportContextForPrompt()  // Get formatted state
```

## 🎯 Benefits for Iterative Development

1. **Accurate Context**: AI gets exact current UI state
2. **Rule Enforcement**: Automatic validation against your standards
3. **Pattern Recognition**: Tracks changes and suggests focus areas
4. **Efficiency**: No manual state description needed
5. **Consistency**: Standardized context format
6. **Proactive**: Detects issues before they become problems

## 📊 Example Use Cases

### Layout Issue Investigation
```javascript
// User reports overlap
window.app.checkRuleCompliance()
// Shows: overlap detected between account-link and object-input

window.app.exportContextForPrompt()
// Provides exact measurements and suggestions
```

### Responsive Design Testing
```javascript
// Resize window to mobile
window.app.captureUIState('mobile_view')

// Check compliance on mobile
window.app.checkRuleCompliance()
// Validates mobile-specific rules
```

### Before/After Comparison
```javascript
// Before changes
window.app.captureUIState('before')

// Make changes...

// After changes
window.app.captureUIState('after')
window.app.exportContextForPrompt()
// Includes change history in context
```

## 🔧 Advanced Features

### Change History Tracking
- Automatically records DOM mutations
- Tracks viewport changes  
- Maintains last 20 changes with timestamps

### Predictive Suggestions
- Analyzes current state for potential issues
- Suggests next actions based on rule compliance
- Prioritizes fixes by impact

### Integration with Existing Tools
- Works with layout debug mode
- Leverages existing .mdc documentation
- Complements error reporting system

This system transforms UI development from reactive debugging to proactive, context-aware iteration. 