import 'package:flutter/material.dart';
import 'package:provider/provider.dart';
import 'payment_service.dart';
import 'camera_service.dart';
import 'analysis_service.dart';
import 'home_screen.dart';
import 'layout_monitor.dart';

void main() {
  // Enable layout monitoring in debug mode
  LayoutMonitor.enable();
  
  runApp(MyApp());
}

class MyApp extends StatelessWidget {
  @override
  Widget build(BuildContext context) {
    return MultiProvider(
      providers: [
        ChangeNotifierProvider(create: (_) => PaymentService()),
        ChangeNotifierProvider(create: (_) => CameraService()),
        ChangeNotifierProvider(create: (_) => AnalysisService()),
      ],
      child: MaterialApp(
        title: 'Molecular Analyzer',
        debugShowCheckedModeBanner: false,
        theme: AppTheme.dark,
        home: HomeScreen(),
        // Add global error handling for layout issues
        builder: (context, child) {
          return _LayoutSafetyWrapper(child: child!);
        },
      ),
    );
  }
}

/// Wrapper to catch and handle layout issues globally
class _LayoutSafetyWrapper extends StatelessWidget {
  final Widget child;
  
  const _LayoutSafetyWrapper({required this.child});
  
  @override
  Widget build(BuildContext context) {
    return MediaQuery(
      data: MediaQuery.of(context),
      child: SafeArea(
        child: LayoutBuilder(
          builder: (context, constraints) {
            // Log screen constraints for debugging
            debugPrint('📐 Screen constraints: ${constraints}');
            
            // Ensure minimum usable space
            if (constraints.maxHeight < 400) {
              return SingleChildScrollView(child: child);
            }
            
            return child;
          },
        ),
      ),
    );
  }
}

/// App theme following ui.mdc guidelines - centralized and safe
class AppTheme {
  static ThemeData get dark => ThemeData.dark().copyWith(
    // Background colors
    scaffoldBackgroundColor: Colors.black,
    canvasColor: Colors.black,
    
    // Card and surface colors
    cardColor: Colors.white.withOpacity(0.08),
    dialogBackgroundColor: Colors.black,
    
    // Primary colors
    primaryColor: Colors.white,
    
    // Text theme
    textTheme: TextTheme(
      bodyLarge: TextStyle(
        color: Colors.white,
        fontSize: 16,
        fontWeight: FontWeight.w400,
        letterSpacing: 0.01,
        height: 1.4,
      ),
      bodyMedium: TextStyle(
        color: Colors.white.withOpacity(0.7),
        fontSize: 14,
        fontWeight: FontWeight.w400,
        letterSpacing: 0.01,
        height: 1.4,
      ),
      titleMedium: TextStyle(
        color: Colors.white,
        fontSize: 18,
        fontWeight: FontWeight.w500,
        letterSpacing: 0.01,
      ),
    ),
    
    // Button themes
    elevatedButtonTheme: ElevatedButtonThemeData(
      style: ElevatedButton.styleFrom(
        backgroundColor: Colors.white.withOpacity(0.08),
        foregroundColor: Colors.white,
        minimumSize: Size(88, 40), // Touch-friendly minimum
        padding: EdgeInsets.symmetric(horizontal: 24, vertical: 16),
        shape: RoundedRectangleBorder(
          borderRadius: BorderRadius.circular(8),
        ),
        textStyle: TextStyle(
          fontSize: 16,
          fontWeight: FontWeight.w500,
          letterSpacing: 0.01,
        ),
      ),
    ),
    
    // Input decoration theme
    inputDecorationTheme: InputDecorationTheme(
      filled: true,
      fillColor: Colors.white.withOpacity(0.08),
      border: OutlineInputBorder(
        borderRadius: BorderRadius.circular(8),
        borderSide: BorderSide.none,
      ),
      hintStyle: TextStyle(
        color: Colors.white.withOpacity(0.5),
        fontSize: 14,
      ),
      contentPadding: EdgeInsets.symmetric(horizontal: 16, vertical: 12),
    ),
    
    // Icon theme
    iconTheme: IconThemeData(
      color: Colors.white.withOpacity(0.7),
      size: 24,
    ),
    
    // Divider theme
    dividerTheme: DividerThemeData(
      color: Colors.white.withOpacity(0.08),
      thickness: 1,
      space: 0,
    ),
  );
} 