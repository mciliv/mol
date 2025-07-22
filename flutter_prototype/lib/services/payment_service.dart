import 'package:flutter/foundation.dart';
import 'package:shared_preferences/shared_preferences.dart';
import '../screens/home_screen.dart'; // For AnalysisMode enum

class PaymentService extends ChangeNotifier {
  bool _showPaymentSection = false;
  bool _isProcessing = false;
  bool _hasPaymentMethod = false;
  String? _userName;
  String? _cardLast4;
  int _usageCount = 0;
  AnalysisMode _selectedMode = AnalysisMode.camera;
  
  // Getters - NO MANUAL DOM QUERIES!
  bool get showPaymentSection => _showPaymentSection;
  bool get isProcessing => _isProcessing;
  bool get hasPaymentMethod => _hasPaymentMethod;
  String? get userName => _userName;
  String? get cardLast4 => _cardLast4;
  int get usageCount => _usageCount;
  AnalysisMode get selectedMode => _selectedMode;
  
  PaymentService() {
    _loadPaymentStatus();
  }
  
  // Load payment status from storage
  Future<void> _loadPaymentStatus() async {
    final prefs = await SharedPreferences.getInstance();
    _hasPaymentMethod = prefs.getBool('hasPaymentMethod') ?? false;
    _userName = prefs.getString('userName');
    _cardLast4 = prefs.getString('cardLast4');
    _usageCount = prefs.getInt('usageCount') ?? 0;
    
    // Show payment section if no payment method
    _showPaymentSection = !_hasPaymentMethod;
    notifyListeners(); // Auto-updates UI - NO MANUAL DOM MANIPULATION!
  }
  
  // Show payment section
  void showPayment() {
    _showPaymentSection = true;
    notifyListeners(); // UI updates automatically!
  }
  
  // Hide payment section
  void hidePayment() {
    _showPaymentSection = false;
    notifyListeners(); // UI updates automatically!
  }
  
  // Setup payment method
  Future<void> setupPayment({String? name}) async {
    _isProcessing = true;
    notifyListeners(); // Button automatically shows loading state
    
    try {
      // Simulate Stripe payment setup
      await Future.delayed(Duration(seconds: 2));
      
      // In real app, this would call Stripe APIs
      // final paymentMethod = await StripeService.createPaymentMethod();
      
      // Store payment info
      final prefs = await SharedPreferences.getInstance();
      await prefs.setBool('hasPaymentMethod', true);
      await prefs.setString('userName', name ?? 'User');
      await prefs.setString('cardLast4', '4242'); // Mock data
      
      // Update state
      _hasPaymentMethod = true;
      _userName = name;
      _cardLast4 = '4242';
      _showPaymentSection = false; // Auto-hide payment section
      
      notifyListeners(); // UI updates automatically!
      
    } catch (error) {
      print('Payment setup failed: $error');
      rethrow;
    } finally {
      _isProcessing = false;
      notifyListeners(); // Button automatically hides loading state
    }
  }
  
  // Change analysis mode
  void setMode(AnalysisMode mode) {
    _selectedMode = mode;
    notifyListeners(); // UI switches automatically!
  }
  
  // Check if analysis is allowed
  bool canAnalyze() {
    return _hasPaymentMethod;
  }
  
  // Increment usage (called after successful analysis)
  Future<void> incrementUsage() async {
    if (!_hasPaymentMethod) return;
    
    _usageCount++;
    notifyListeners(); // UI shows updated count automatically
    
    // Save to storage
    final prefs = await SharedPreferences.getInstance();
    await prefs.setInt('usageCount', _usageCount);
    
    // In real app, also send to backend
    // await ApiService.incrementUsage();
  }
  
  // Setup developer account (for testing)
  Future<void> setupDeveloperAccount() async {
    final prefs = await SharedPreferences.getInstance();
    await prefs.setBool('hasPaymentMethod', true);
    await prefs.setString('userName', 'Developer');
    await prefs.setString('cardLast4', '0000');
    
    _hasPaymentMethod = true;
    _userName = 'Developer';
    _cardLast4 = '0000';
    _showPaymentSection = false;
    
    notifyListeners(); // UI updates automatically!
  }
  
  // Reset payment (for testing)
  Future<void> resetPayment() async {
    final prefs = await SharedPreferences.getInstance();
    await prefs.clear();
    
    _hasPaymentMethod = false;
    _userName = null;
    _cardLast4 = null;
    _usageCount = 0;
    _showPaymentSection = true;
    
    notifyListeners(); // UI resets automatically!
  }
} 