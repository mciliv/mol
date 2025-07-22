import 'package:flutter/foundation.dart';
import 'dart:io';

class AnalysisService extends ChangeNotifier {
  bool _isAnalyzing = false;
  List<MolecularResult> _results = [];
  String? _errorMessage;
  
  // Getters
  bool get isAnalyzing => _isAnalyzing;
  List<MolecularResult> get results => _results;
  String? get errorMessage => _errorMessage;
  
  // Analyze text input
  Future<void> analyzeText(String input) async {
    if (input.trim().isEmpty) return;
    
    _isAnalyzing = true;
    _errorMessage = null;
    notifyListeners();
    
    try {
      // Simulate molecular analysis
      await Future.delayed(Duration(seconds: 2));
      
      // Mock result - in real app, this would call your backend
      final result = MolecularResult(
        name: input,
        formula: 'C8H10N4O2', // Caffeine as example
        description: 'Analysis result for $input',
        timestamp: DateTime.now(),
      );
      
      _results.insert(0, result);
      notifyListeners();
      
    } catch (error) {
      _errorMessage = error.toString();
      notifyListeners();
    } finally {
      _isAnalyzing = false;
      notifyListeners();
    }
  }
  
  // Analyze image file
  Future<void> analyzeImage(File imageFile) async {
    _isAnalyzing = true;
    _errorMessage = null;
    notifyListeners();
    
    try {
      // Simulate image analysis
      await Future.delayed(Duration(seconds: 3));
      
      // Mock result - in real app, this would use OpenAI Vision API
      final result = MolecularResult(
        name: 'Image Analysis Result',
        formula: 'C6H12O6', // Glucose as example
        description: 'Molecular analysis from captured image',
        timestamp: DateTime.now(),
        imagePath: imageFile.path,
      );
      
      _results.insert(0, result);
      notifyListeners();
      
    } catch (error) {
      _errorMessage = error.toString();
      notifyListeners();
    } finally {
      _isAnalyzing = false;
      notifyListeners();
    }
  }
  
  // Analyze image from URL
  Future<void> analyzeImageUrl(String url) async {
    if (url.trim().isEmpty) return;
    
    _isAnalyzing = true;
    _errorMessage = null;
    notifyListeners();
    
    try {
      // Simulate URL image analysis
      await Future.delayed(Duration(seconds: 2));
      
      final result = MolecularResult(
        name: 'URL Image Analysis',
        formula: 'H2O', // Water as example
        description: 'Analysis from image URL: $url',
        timestamp: DateTime.now(),
        imageUrl: url,
      );
      
      _results.insert(0, result);
      notifyListeners();
      
    } catch (error) {
      _errorMessage = error.toString();
      notifyListeners();
    } finally {
      _isAnalyzing = false;
      notifyListeners();
    }
  }
  
  // Clear results
  void clearResults() {
    _results.clear();
    notifyListeners();
  }
  
  // Remove specific result
  void removeResult(int index) {
    if (index >= 0 && index < _results.length) {
      _results.removeAt(index);
      notifyListeners();
    }
  }
}

class MolecularResult {
  final String name;
  final String formula;
  final String description;
  final DateTime timestamp;
  final String? imagePath;
  final String? imageUrl;
  
  MolecularResult({
    required this.name,
    required this.formula,
    required this.description,
    required this.timestamp,
    this.imagePath,
    this.imageUrl,
  });
} 