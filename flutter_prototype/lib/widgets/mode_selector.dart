import 'package:flutter/material.dart';
import 'package:provider/provider.dart';
import '../services/payment_service.dart';
import '../screens/home_screen.dart'; // For AnalysisMode enum

class ModeSelector extends StatelessWidget {
  @override
  Widget build(BuildContext context) {
    return Consumer<PaymentService>(
      builder: (context, paymentService, child) {
        return Container(
          padding: EdgeInsets.all(8),
          decoration: BoxDecoration(
            color: Colors.white.withOpacity(0.08), // Interactive element background
            borderRadius: BorderRadius.circular(8),
          ),
          child: Row(
            children: [
              Expanded(
                child: _buildModeTab(
                  context,
                  paymentService,
                  AnalysisMode.camera,
                  Icons.camera_alt,
                  'Camera',
                ),
              ),
              SizedBox(width: 8),
              Expanded(
                child: _buildModeTab(
                  context,
                  paymentService,
                  AnalysisMode.photo,
                  Icons.photo,
                  'Photo',
                ),
              ),
            ],
          ),
        );
      },
    );
  }
  
  Widget _buildModeTab(
    BuildContext context,
    PaymentService paymentService,
    AnalysisMode mode,
    IconData icon,
    String label,
  ) {
    final isSelected = paymentService.selectedMode == mode;
    
    return GestureDetector(
      onTap: () => paymentService.setMode(mode),
      child: Container(
        padding: EdgeInsets.symmetric(vertical: 12, horizontal: 16),
        decoration: BoxDecoration(
          color: isSelected 
            ? Colors.white.withOpacity(0.08) 
            : Colors.transparent,
          borderRadius: BorderRadius.circular(8),
        ),
        child: Row(
          mainAxisAlignment: MainAxisAlignment.center,
          children: [
            Icon(
              icon,
              size: 16,
              color: isSelected 
                ? Colors.white 
                : Colors.white.withOpacity(0.7),
            ),
            SizedBox(width: 8),
            Text(
              label,
              style: TextStyle(
                color: isSelected 
                  ? Colors.white 
                  : Colors.white.withOpacity(0.7),
                fontWeight: FontWeight.w500,
                fontSize: 14,
              ),
            ),
          ],
        ),
      ),
    );
  }
} 