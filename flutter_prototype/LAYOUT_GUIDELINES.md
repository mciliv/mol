# Flutter Layout Guidelines 📐

## 🚨 **Critical Layout Rules - NEVER IGNORE**

### 1. **Column/Row Overflow Prevention**
```dart
// ❌ WRONG - Can cause overflow
Column(
  children: [
    Text('Long text that might overflow'),
    Container(height: 200, child: SomeWidget()),
    // More widgets...
  ],
)

// ✅ CORRECT - Prevents overflow
Column(
  children: [
    Text('Long text that might overflow').flexible(),
    Container(height: 200, child: SomeWidget()).expanded(),
    // Fixed-size widgets last
  ],
)

// ✅ EVEN BETTER - Use our utilities
LayoutMonitor.safeFlex(
  children: [
    Text('Safe text'),
    SomeWidget(),
  ],
  allowOverflow: true, // Adds scrolling if needed
)
```

### 2. **Always Use MainAxisSize.min for Nested Columns**
```dart
// ❌ WRONG - Takes all available space
Column(
  children: [Icon(), Text()],
)

// ✅ CORRECT - Only takes needed space
Column(
  mainAxisSize: MainAxisSize.min, // KEY ADDITION
  children: [Icon(), Text()],
)
```

### 3. **Fixed Heights for Known Content**
```dart
// ❌ WRONG - Unbounded height
Container(
  child: Column(children: [/* dynamic content */]),
)

// ✅ CORRECT - Bounded height
Container(
  height: 200, // or use ConstrainedBox
  child: Column(children: [/* dynamic content */]),
)
```

## 🛠️ **Required Tools Usage**

### Use Layout Extensions
```dart
import '../utils/layout_monitor.dart';

// Prevent overflow automatically
Widget myWidget = SomeComplexWidget()
  .preventOverflow(maxHeight: 300)
  .scrollable(); // Add scrolling if needed

// Safe flexible layouts
Widget safeLayout = Column(
  children: [
    TopWidget().flexible(),
    MiddleWidget().expanded(flex: 2),
    BottomWidget().flexible(),
  ],
);
```

### Use Safe Containers
```dart
// Instead of regular Container
LayoutMonitor.safeContainer(
  maxHeight: 400, // Automatic constraint
  padding: EdgeInsets.all(16),
  child: YourContent(),
)
```

## 📱 **Screen Size Guidelines**

### Responsive Design Rules
```dart
// Use LayoutBuilder for different screen sizes
LayoutBuilder(
  builder: (context, constraints) {
    if (constraints.maxWidth < 600) {
      return MobileLayout();
    } else {
      return DesktopLayout();
    }
  },
)

// Or use MediaQuery for breakpoints
final screenWidth = MediaQuery.of(context).size.width;
final isMobile = screenWidth < 600;
```

### Fixed Dimensions Guidelines
- **Minimum touch target**: 40x40 pixels
- **Maximum modal height**: 80% of screen height
- **Fixed widget heights**: Use multiples of 20 (40, 60, 80, 100, 120...)
- **Container padding**: Use 8, 16, 20, 24, 32 for consistency

## 🔧 **Widget-Specific Rules**

### TextFields
```dart
// ✅ Always constrain in forms
ConstrainedBox(
  constraints: BoxConstraints(maxWidth: 400),
  child: TextField(
    decoration: InputDecoration(
      isDense: true, // Reduces height
      contentPadding: EdgeInsets.symmetric(
        horizontal: 16, 
        vertical: 12,
      ),
    ),
  ),
)
```

### Lists and Scrolling
```dart
// ✅ Always bound ListView in Column
Column(
  children: [
    SomeHeader(),
    Expanded( // KEY: Gives remaining space to list
      child: ListView.builder(
        itemCount: items.length,
        itemBuilder: (context, index) => ItemWidget(items[index]),
      ),
    ),
  ],
)
```

### Images and Media
```dart
// ✅ Always constrain images
Image.asset(
  'path/to/image.png',
  width: double.infinity,
  height: 200, // Fixed height
  fit: BoxFit.cover, // Prevents distortion
)
```

## 🧪 **Testing Requirements**

### Widget Tests
```dart
// Test different screen sizes
testWidgets('Layout works on mobile', (tester) async {
  // Set mobile screen size
  await tester.binding.setSurfaceSize(Size(375, 667));
  await tester.pumpWidget(MyApp());
  
  // Verify no overflow
  expect(tester.takeException(), isNull);
});

testWidgets('Layout works on tablet', (tester) async {
  // Set tablet screen size
  await tester.binding.setSurfaceSize(Size(768, 1024));
  await tester.pumpWidget(MyApp());
  
  // Verify no overflow
  expect(tester.takeException(), isNull);
});
```

### Golden Tests for Layout
```dart
// Create golden tests for critical layouts
testWidgets('Payment section golden test', (tester) async {
  await tester.pumpWidget(PaymentSection());
  await expectLater(
    find.byType(PaymentSection),
    matchesGoldenFile('payment_section.png'),
  );
});
```

## 🚦 **Development Workflow**

### Pre-Commit Checklist
- [ ] Run `flutter analyze` - no layout warnings
- [ ] Test on mobile size (375x667)
- [ ] Test on tablet size (768x1024)
- [ ] Test on desktop size (1200x800)
- [ ] Check console for overflow errors
- [ ] Verify all touch targets ≥ 40px

### Debug Commands
```bash
# Check for layout issues
flutter run --debug

# Enable layout debugging
flutter run --debug --enable-software-rendering

# Performance debugging
flutter run --profile
```

## 🎯 **Common Fixes**

### Overflow Errors
1. **Immediate fix**: Wrap in `SingleChildScrollView`
2. **Better fix**: Use `Flexible/Expanded` widgets
3. **Best fix**: Use our `LayoutMonitor.safeFlex`

### Performance Issues
1. **Use ListView for many items** (>20)
2. **Const constructors** for static widgets
3. **RepaintBoundary** for complex widgets

### Accessibility
1. **Semantic labels** for all interactive elements
2. **Minimum contrast ratios** (4.5:1 for normal text)
3. **Touch targets ≥ 40px**

## 🔄 **Continuous Monitoring**

The `LayoutMonitor` utility will:
- ✅ Catch overflow errors automatically
- ✅ Provide helpful fix suggestions
- ✅ Log screen constraints for debugging
- ✅ Auto-wrap problematic widgets

**Remember**: Layout issues caught early are 10x easier to fix than after deployment! 🎯 