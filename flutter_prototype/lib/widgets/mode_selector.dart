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
                  _buildCameraIcon(),
                  'Camera',
                ),
              ),
              SizedBox(width: 8),
              Expanded(
                child: _buildModeTab(
                  context,
                  paymentService,
                  AnalysisMode.photo,
                  _buildPhotoIcon(),
                  'Photo',
                ),
              ),
            ],
          ),
        );
      },
    );
  }

  // Original camera SVG from HTML: camera body with lens
  Widget _buildCameraIcon() {
    return CustomPaint(
      size: Size(16, 16),
      painter: CameraIconPainter(),
    );
  }

  // Original photo SVG from HTML: landscape/mountain scene
  Widget _buildPhotoIcon() {
    return CustomPaint(
      size: Size(16, 16),
      painter: PhotoIconPainter(),
    );
  }
  
  Widget _buildModeTab(
    BuildContext context,
    PaymentService paymentService,
    AnalysisMode mode,
    Widget icon,
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
            Theme(
              data: Theme.of(context).copyWith(
                iconTheme: IconThemeData(
                  color: isSelected 
                    ? Colors.white 
                    : Colors.white.withOpacity(0.7),
                ),
              ),
              child: icon,
            ),
            SizedBox(width: 8),
            Text(
              label,
              style: TextStyle(
                color: isSelected 
                  ? Colors.white 
                  : Colors.white.withOpacity(0.7),
                fontSize: 14,
                fontWeight: FontWeight.w500,
              ),
            ),
          ],
        ),
      ),
    );
  }
}

// Custom painter for camera icon matching original SVG
class CameraIconPainter extends CustomPainter {
  @override
  void paint(Canvas canvas, Size size) {
    final paint = Paint()
      ..color = Colors.white
      ..style = PaintingStyle.stroke
      ..strokeWidth = 2.0;

    final path = Path();
    
    // Scale factor to fit 16x16 canvas from 24x24 viewBox
    final scale = size.width / 24.0;
    
    // Camera body: M14.5 4h-5L7 6H4a2 2 0 0 0-2 2v9a2 2 0 0 0 2 2h16a2 2 0 0 0 2-2V8a2 2 0 0 0-2-2h-3l-2.5-2z
    path.moveTo(14.5 * scale, 4 * scale);
    path.lineTo(9.5 * scale, 4 * scale);
    path.lineTo(7 * scale, 6 * scale);
    path.lineTo(4 * scale, 6 * scale);
    path.addRRect(RRect.fromLTRBR(
      2 * scale, 6 * scale, 22 * scale, 17 * scale,
      Radius.circular(2 * scale)
    ));
    path.moveTo(19 * scale, 6 * scale);
    path.lineTo(16.5 * scale, 4 * scale);
    
    canvas.drawPath(path, paint);
    
    // Camera lens: circle cx="12" cy="13" r="3"
    canvas.drawCircle(
      Offset(12 * scale, 13 * scale),
      3 * scale,
      paint,
    );
  }

  @override
  bool shouldRepaint(CustomPainter oldDelegate) => false;
}

// Custom painter for photo icon matching original SVG
class PhotoIconPainter extends CustomPainter {
  @override
  void paint(Canvas canvas, Size size) {
    final paint = Paint()
      ..color = Colors.white
      ..style = PaintingStyle.stroke
      ..strokeWidth = 2.0;

    final scale = size.width / 24.0;
    
    // Photo frame: rect x="3" y="3" width="18" height="18" rx="2" ry="2"
    canvas.drawRRect(
      RRect.fromLTRBR(
        3 * scale, 3 * scale, 21 * scale, 21 * scale,
        Radius.circular(2 * scale)
      ),
      paint,
    );
    
    // Sun/moon: circle cx="9" cy="9" r="2"
    canvas.drawCircle(
      Offset(9 * scale, 9 * scale),
      2 * scale,
      paint,
    );
    
    // Mountain/landscape path: M21 15l-3.086-3.086a2 2 0 0 0-2.828 0L6 21
    final mountainPath = Path();
    mountainPath.moveTo(21 * scale, 15 * scale);
    mountainPath.lineTo(17.914 * scale, 11.914 * scale);
    mountainPath.quadraticBezierTo(
      16.5 * scale, 10.5 * scale,
      15.086 * scale, 11.914 * scale
    );
    mountainPath.lineTo(6 * scale, 21 * scale);
    
    canvas.drawPath(mountainPath, paint);
  }

  @override
  bool shouldRepaint(CustomPainter oldDelegate) => false;
} 