import 'package:flutter/material.dart';
import 'package:provider/provider.dart';
import '../widgets/payment_section.dart';
import '../widgets/input_bar.dart';
import '../widgets/mode_selector.dart';
import '../widgets/camera_view.dart';
import '../widgets/photo_picker.dart';
import '../widgets/results_view.dart';
import '../services/payment_service.dart';

class HomeScreen extends StatelessWidget {
  @override
  Widget build(BuildContext context) {
    return Scaffold(
      backgroundColor: Colors.black,
      body: SafeArea(
        child: Consumer<PaymentService>(
          builder: (context, paymentService, child) {
            return Column(
              children: [
                // Payment section (shows when needed) - NO MODAL COMPLEXITY!
                if (paymentService.showPaymentSection) ...[
                  PaymentSection(),
                  Container(
                    height: 1,
                    color: Colors.white.withOpacity(0.08), // Simple separator
                  ),
                ],
                
                // Main app content
                Expanded(
                  child: Padding(
                    padding: EdgeInsets.all(20),
                    child: Column(
                      children: [
                        InputBar(), // Text input + account button
                        SizedBox(height: 20),
                        ModeSelector(), // Camera/Photo tabs
                        SizedBox(height: 20),
                        Expanded(child: ContentArea()), // Camera/Photo/Results
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

class ContentArea extends StatelessWidget {
  @override
  Widget build(BuildContext context) {
    return Consumer<PaymentService>(
      builder: (context, paymentService, child) {
        // Simple state-based UI switching
        switch (paymentService.selectedMode) {
          case AnalysisMode.camera:
            return CameraView();
          case AnalysisMode.photo:
            return PhotoPicker();
          default:
            return ResultsView();
        }
      },
    );
  }
}

enum AnalysisMode { camera, photo, results } 