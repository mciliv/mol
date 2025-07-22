import 'package:flutter/material.dart';
import 'package:flutter/rendering.dart';

/// Layout monitoring utility to prevent overflow and constraint issues
class LayoutMonitor {
  static bool _isEnabled = false;
  
  /// Enable layout monitoring in debug mode
  static void enable() {
    assert(() {
      _isEnabled = true;
      // Override Flutter's error handler for layout issues
      FlutterError.onError = (FlutterErrorDetails details) {
        if (details.exception.toString().contains('RenderFlex overflowed')) {
          _handleLayoutOverflow(details);
        }
        // Still call the original handler
        FlutterError.presentError(details);
      };
      return true;
    }());
  }
  
  /// Handle layout overflow errors with helpful debugging
  static void _handleLayoutOverflow(FlutterErrorDetails details) {
    print('\n🚨 LAYOUT OVERFLOW DETECTED:');
    print('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
    print('❌ ${details.exception}');
    print('📍 Widget: ${details.context}');
    print('\n💡 SOLUTIONS:');
    print('• Use Expanded/Flexible widgets');
    print('• Add SingleChildScrollView');
    print('• Set fixed heights with ConstrainedBox');
    print('• Use mainAxisSize: MainAxisSize.min');
    print('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');
  }
  
  /// Wrap widgets to prevent common layout issues
  static Widget safeFlex({
    required List<Widget> children,
    Axis direction = Axis.vertical,
    MainAxisAlignment mainAxisAlignment = MainAxisAlignment.start,
    CrossAxisAlignment crossAxisAlignment = CrossAxisAlignment.center,
    bool allowOverflow = false,
  }) {
    if (allowOverflow) {
      return SingleChildScrollView(
        scrollDirection: direction,
        child: Flex(
          direction: direction,
          mainAxisAlignment: mainAxisAlignment,
          crossAxisAlignment: crossAxisAlignment,
          children: children,
        ),
      );
    }
    
    return Flex(
      direction: direction,
      mainAxisAlignment: mainAxisAlignment,
      crossAxisAlignment: crossAxisAlignment,
      children: children.map((child) {
        // Auto-wrap in Flexible to prevent overflow
        if (child is Expanded || child is Flexible) {
          return child;
        }
        return Flexible(child: child);
      }).toList(),
    );
  }
  
  /// Safe container with automatic constraints
  static Widget safeContainer({
    required Widget child,
    double? width,
    double? height,
    double? maxWidth,
    double? maxHeight,
    EdgeInsetsGeometry? padding,
    EdgeInsetsGeometry? margin,
    Decoration? decoration,
  }) {
    return ConstrainedBox(
      constraints: BoxConstraints(
        maxWidth: maxWidth ?? double.infinity,
        maxHeight: maxHeight ?? double.infinity,
      ),
      child: Container(
        width: width,
        height: height,
        padding: padding,
        margin: margin,
        decoration: decoration,
        child: child,
      ),
    );
  }
}

/// Extension to add layout safety to any widget
extension LayoutSafety on Widget {
  /// Wrap widget to prevent overflow
  Widget preventOverflow({double? maxHeight, double? maxWidth}) {
    return ConstrainedBox(
      constraints: BoxConstraints(
        maxWidth: maxWidth ?? double.infinity,
        maxHeight: maxHeight ?? double.infinity,
      ),
      child: this,
    );
  }
  
  /// Wrap in flexible to play nice with Flex widgets
  Widget flexible({int flex = 1}) {
    return Flexible(flex: flex, child: this);
  }
  
  /// Wrap in expanded for flex layouts
  Widget expanded({int flex = 1}) {
    return Expanded(flex: flex, child: this);
  }
  
  /// Add safe scrolling if content might overflow
  Widget scrollable({Axis scrollDirection = Axis.vertical}) {
    return SingleChildScrollView(
      scrollDirection: scrollDirection,
      child: this,
    );
  }
} 