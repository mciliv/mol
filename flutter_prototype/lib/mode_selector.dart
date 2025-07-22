import 'package:flutter/material.dart';
import 'package:provider/provider.dart';
import 'package:flutter_svg/flutter_svg.dart';
import 'payment_service.dart';

class ModeSelector extends StatelessWidget {
  @override
  Widget build(BuildContext context) {
    return Consumer<PaymentService>(
      builder: (context, paymentService, child) {
        return Row(
          children: [
            Expanded(
              child: _buildModeTab(
                context,
                paymentService,
                AnalysisMode.camera,
                'Camera',
                'assets/lens.svg',
              ),
            ),
            SizedBox(width: 12),
            Expanded(
              child: _buildModeTab(
                context,
                paymentService,
                AnalysisMode.photo,
                'Photo',
                'assets/photo-upload.svg',
              ),
            ),
          ],
        );
      },
    );
  }
  
  Widget _buildModeTab(
    BuildContext context,
    PaymentService paymentService,
    AnalysisMode mode,
    String label,
    String svgAsset,
  ) {
    final isSelected = paymentService.selectedMode == mode;
    
    return GestureDetector(
      onTap: () => paymentService.setSelectedMode(mode),
      child: Container(
        padding: EdgeInsets.symmetric(horizontal: 16, vertical: 12),
        decoration: BoxDecoration(
          color: isSelected 
            ? Colors.white.withOpacity(0.12) 
            : Colors.white.withOpacity(0.08),
          borderRadius: BorderRadius.circular(8),
        ),
        child: Row(
          mainAxisAlignment: MainAxisAlignment.center,
          children: [
            SvgPicture.asset(
              svgAsset,
              width: 16,
              height: 16,
              colorFilter: ColorFilter.mode(
                isSelected ? Colors.white : Colors.white.withOpacity(0.7),
                BlendMode.srcIn,
              ),
            ),
            SizedBox(width: 8),
            Text(
              label,
              style: TextStyle(
                color: isSelected ? Colors.white : Colors.white.withOpacity(0.7),
                fontWeight: isSelected ? FontWeight.w600 : FontWeight.w500,
                fontSize: 14,
              ),
            ),
          ],
        ),
      ),
    );
  }
} 