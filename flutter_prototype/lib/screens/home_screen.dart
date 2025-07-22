import 'package:flutter/material.dart';
import 'package:provider/provider.dart';
import '../services/payment_service.dart';
import '../widgets/payment_section.dart';
import '../widgets/input_bar.dart';
import '../widgets/mode_selector.dart';
import '../widgets/photo_picker.dart';
import '../widgets/results_view.dart';

class HomeScreen extends StatelessWidget {
  @override
  Widget build(BuildContext context) {
    return Scaffold(
      backgroundColor: Colors.black, // Following ui.mdc
      body: SafeArea(
        child: Consumer<PaymentService>(
          builder: (context, paymentService, child) {
            return Column(
              children: [
                // Payment Section (conditionally shown)
                if (paymentService.showPaymentSection) PaymentSection(),
                
                // Separator Line (conditionally shown)
                if (paymentService.showPaymentSection)
                  Container(
                    height: 1,
                    color: Colors.white.withOpacity(0.08),
                  ),
                
                // Main App Content - CLEAN & FLOWING (not boxy)
                Expanded(
                  child: SingleChildScrollView(
                    padding: EdgeInsets.all(20),
                    child: Column(
                      children: [
                        // Input Bar
                        InputBar(),
                        SizedBox(height: 20),
                        
                        // Mode Selector
                        ModeSelector(),
                        SizedBox(height: 20),
                        
                        // Mode-specific content - FLOWING (no fixed heights)
                        _buildModeContent(paymentService),
                        
                        SizedBox(height: 20),
                        
                        // Results View - MINIMAL CONSTRAINTS
                        SizedBox(
                          height: 300, // Sufficient space, not visually constraining
                          child: ResultsView(),
                        ),
                      ],
                    ),
                  ),
                ),
              ],
            );
          },
        ),
      ),
    );
  }
  
  // Clean mode content - no visual boxes
  Widget _buildModeContent(PaymentService paymentService) {
    switch (paymentService.selectedMode) {
      case AnalysisMode.photo:
        return PhotoPicker();
      case AnalysisMode.camera:
        return SimpleCameraView();
      default:
        return SizedBox.shrink();
    }
  }
}

// Clean Camera View - no boxy container
class SimpleCameraView extends StatelessWidget {
  @override
  Widget build(BuildContext context) {
    return Padding(
      padding: EdgeInsets.symmetric(vertical: 40),
      child: Row(
        mainAxisAlignment: MainAxisAlignment.center,
        children: [
          Icon(
            Icons.camera_alt,
            size: 32,
            color: Colors.white.withOpacity(0.3),
          ),
          SizedBox(width: 16),
          Column(
            mainAxisSize: MainAxisSize.min, // KEY: Prevents overflow
            children: [
              Text(
                'Camera View',
                style: TextStyle(
                  color: Colors.white.withOpacity(0.7),
                  fontSize: 16,
                ),
              ),
              Text(
                'Coming soon',
                style: TextStyle(
                  color: Colors.white.withOpacity(0.5),
                  fontSize: 12,
                ),
              ),
            ],
          ),
        ],
      ),
    );
  }
}

enum AnalysisMode { camera, photo, results } 