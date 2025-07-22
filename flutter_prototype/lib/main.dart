import 'package:flutter/material.dart';
import 'package:provider/provider.dart';
import 'screens/home_screen.dart';
import 'services/payment_service.dart';
import 'services/analysis_service.dart';
import 'services/camera_service.dart';

void main() {
  runApp(MolecularAnalyzerApp());
}

class MolecularAnalyzerApp extends StatelessWidget {
  @override
  Widget build(BuildContext context) {
    return MultiProvider(
      providers: [
        ChangeNotifierProvider(create: (_) => PaymentService()),
        ChangeNotifierProvider(create: (_) => AnalysisService()),
        ChangeNotifierProvider(create: (_) => CameraService()),
      ],
      child: MaterialApp(
        title: 'Molecular Analyzer',
        debugShowCheckedModeBanner: false,
        theme: AppTheme.darkTheme,
        home: HomeScreen(),
      ),
    );
  }
}

class AppTheme {
  static final ThemeData darkTheme = ThemeData.dark().copyWith(
    // Following ui.mdc guidelines - "stupid simple" dark theme
    scaffoldBackgroundColor: Colors.black,
    cardColor: Colors.white.withOpacity(0.08),
    
    // Clean button styling
    elevatedButtonTheme: ElevatedButtonThemeData(
      style: ElevatedButton.styleFrom(
        backgroundColor: Colors.white.withOpacity(0.08),
        foregroundColor: Colors.white,
        minimumSize: Size(double.infinity, 40), // Touch-friendly
        shape: RoundedRectangleBorder(borderRadius: BorderRadius.circular(8)),
        elevation: 0, // No extraneous shadows
      ),
    ),
    
    // Simple input styling
    inputDecorationTheme: InputDecorationTheme(
      filled: true,
      fillColor: Colors.white.withOpacity(0.08),
      border: OutlineInputBorder(
        borderRadius: BorderRadius.circular(8),
        borderSide: BorderSide.none, // No extraneous lines
      ),
      hintStyle: TextStyle(color: Colors.white.withOpacity(0.7)),
    ),
    
    // Clean text theme
    textTheme: TextTheme(
      bodyLarge: TextStyle(
        fontFamily: '-apple-system', // Following ui.mdc typography
        fontWeight: FontWeight.w400,
        letterSpacing: 0.01,
        height: 1.4,
        color: Colors.white,
      ),
      headlineSmall: TextStyle(
        fontWeight: FontWeight.w500, // Label weight
        color: Colors.white,
      ),
    ),
  );
} 