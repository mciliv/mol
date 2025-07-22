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
                
                // Main App Content - FLATTENED STRUCTURE
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
                        
                        // Mode-specific content - SIMPLIFIED
                        _buildModeContent(paymentService),
                        
                        SizedBox(height: 20),
                        
                        // Results View - CONSTRAINED HEIGHT
                        ConstrainedBox(
                          constraints: BoxConstraints(
                            minHeight: 200,
                            maxHeight: 400,
                          ),
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
  
  // Flattened mode content - no complex nesting
  Widget _buildModeContent(PaymentService paymentService) {
    switch (paymentService.selectedMode) {
      case AnalysisMode.photo:
        return Container(
          height: 200, // Fixed height prevents overflow
          child: PhotoPicker(),
        );
      case AnalysisMode.camera:
        return Container(
          height: 200, // Fixed height prevents overflow
          child: SimpleCameraView(),
        );
      default:
        return SizedBox.shrink();
    }
  }
}

// Simplified Camera View - no nested columns
class SimpleCameraView extends StatelessWidget {
  @override
  Widget build(BuildContext context) {
    return Container(
      width: double.infinity,
      decoration: BoxDecoration(
        color: Colors.white.withOpacity(0.08),
        borderRadius: BorderRadius.circular(8),
      ),
      child: Center(
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
      ),
    );
  }
}

enum AnalysisMode { camera, photo, results } 