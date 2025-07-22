import 'package:flutter/material.dart';
import 'package:provider/provider.dart';
import '../services/payment_service.dart';
import '../services/camera_service.dart';
import '../services/analysis_service.dart';
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
                
                // Main App Content
                Expanded(
                  child: Padding(
                    padding: EdgeInsets.all(20),
                    child: Column(
                      children: [
                        // Input Bar
                        InputBar(),
                        SizedBox(height: 20),
                        
                        // Mode Selector
                        ModeSelector(),
                        SizedBox(height: 20),
                        
                        // Content Area - Now includes results
                        Expanded(
                          child: Column(
                            children: [
                              // Mode-specific content (camera/photo)
                              if (paymentService.selectedMode == AnalysisMode.photo)
                                Expanded(flex: 1, child: PhotoPicker()),
                              
                              if (paymentService.selectedMode == AnalysisMode.camera)
                                Expanded(flex: 1, child: CameraView()),
                              
                              SizedBox(height: 20),
                              
                              // Results View - Always shown
                              Expanded(
                                flex: 2, // Give more space to results
                                child: ResultsView(),
                              ),
                            ],
                          ),
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
}

// Camera placeholder widget
class CameraView extends StatelessWidget {
  @override
  Widget build(BuildContext context) {
    return Container(
      decoration: BoxDecoration(
        color: Colors.white.withOpacity(0.08),
        borderRadius: BorderRadius.circular(8),
      ),
      child: Center(
        child: Column(
          mainAxisAlignment: MainAxisAlignment.center,
          children: [
            Icon(
              Icons.camera_alt,
              size: 64,
              color: Colors.white.withOpacity(0.3),
            ),
            SizedBox(height: 16),
            Text(
              'Camera View',
              style: TextStyle(
                color: Colors.white.withOpacity(0.7),
                fontSize: 16,
              ),
            ),
            SizedBox(height: 8),
            Text(
              'Camera integration coming soon',
              style: TextStyle(
                color: Colors.white.withOpacity(0.5),
                fontSize: 12,
              ),
            ),
          ],
        ),
      ),
    );
  }
}

enum AnalysisMode { camera, photo, results } 